# This is code Op_en3104/pollutant_bayes on page [[EU-kalat]]
# The code is also available at https://github.com/jtuomist/pfas/blob/main/conc_pcddf_preprocess.R

library(OpasnetUtils)
library(reshape2)
library(rjags) # JAGS
library(ggplot2)
library(MASS) # mvrnorm
library(car) # scatterplotMatrix

#size <- Ovariable("size", ddata="Op_en7748", subset="Size distribution of fish species")
#time <- Ovariable("time", data = data.frame(Result=2015))
#conc_pcddf <- EvalOutput(conc_pcddf,verbose=TRUE)
#View(conc_pcddf@output)

objects.latest("Op_en3104", code_name = "preprocess") # [[EU-kalat]] eu, eu2, euRatio, indices

eu2 <- EvalOutput(eu2)

# Hierarchical Bayes model.

# PCDD/F concentrations in fish.
# It uses the TEQ sum of PCDD/F (PCDDF) as the total concentration
# of dioxin and PCB respectively for PCB in fish.
# PCDDF depends on size of fish, fish species, catchment time, and catchment area,
# but we omit catchment area. In addition, we assume that size of fish has
# zero effect for other fish than Baltic herring.
# Catchment year affects all species similarly. 

eu2 <- eu2[!eu2$Compound %in% c("MPhT","DOT","BDE138"),] # No values > 0

eu3 <- eu2[eu2$Matrix == "Muscle" , ]@output
eu3 <- reshape( 
  eu3, 
  v.names = "eu2Result", 
  idvar = c("THLcode", "Fish"),
  timevar = "Compound", 
  drop = c("Matrix","eu2Source"), 
  direction = "wide"
)
colnames(eu3) <- gsub("eu2Result\\.","",colnames(eu3))

conl_nd <- c("PFOA","PFOS","DBT","MBT","TBT","DPhT","TPhT")
eu4 <- eu3[rowSums(is.na(eu3[conl_nd]))<7 , c(1:5,match(conl_nd,colnames(eu3)))]
fisl_nd <- as.character(unique(eu4$Fish))
conc_nd <- data.matrix(eu4[6:ncol(eu4)])

eu3 <- eu3[!is.na(eu3$PCDDF) , !colnames(eu3) %in% conl_nd]
conl <- colnames(eu3)[-(1:5)]
fisl <- as.character(unique(eu3$Fish))

oprint(head(eu3))

C <- length(conl)
Fi <- length(fisl)
N <- 1000
conl
fisl

conc <- eu3[6:ncol(eu3)]

# Find the level of quantification for dinterval function
LOQ <- unlist(lapply(eu3[6:ncol(eu3)], FUN = function(x) min(x[x!=0], na.rm=TRUE))) 
# With TEQ, there are no zeros. So this is useful only if there are congener-specific results.
#names(LOQ) <- conl


conc <- sapply(
  1:length(LOQ),
  FUN = function(x) ifelse(is.na(conc[,x]) | conc[,x]==0, 0.5*LOQ[x], conc[,x])
)
conc <- data.matrix(conc)

# It assumes that all fish groups have the same Omega but mu varies.

mod <- textConnection(
  "
  model{
    for(i in 1:S) { # S = fish sample
      #        below.LOQ[i,j] ~ dinterval(-conc[i,j], -LOQ[j])
      conc[i,1:C] ~ dmnorm(muind[i,], Omega[fis[i],,])
      muind[i,1:C] <- mu[fis[i],1:C] #+ lenp[fis[i]]*length[i] + timep*year[i]
    }
    for(i in 1:S_nd) {
      for(j in 1:C_nd) {
        conc[i,j] ~ dnorm(muind_nd[i,], tau_nd[fis_nd[i],])
        muind[i,j] <- mu[fis_nd[i],j] #+ lenp[fis[i]]*length[i] + timep*year[i]
      }
    }
    
    # Priors for parameters
    # Time trend. Assumed a known constant because at the moment there is little time variation in data.
    # https://www.evira.fi/elintarvikkeet/ajankohtaista/2018/itameren-silakoissa-yha-vahemman-ymparistomyrkkyja---paastojen-rajoitukset-vaikuttavat/
    # PCDDF/PCB-concentations 2001: 9 pg/g fw, 2016: 3.5 pg/g fw. (3.5/9)^(1/15)-1=-0.06102282
  #  timep ~ dnorm(-0.0610, 10000) 
  #  lenp[1] ~ dnorm(0.01,0.01) # length parameter for herring
  #  lenp[2] ~ dnorm(0,10000) # length parameter for salmon: assumed zero
    
    for(i in 1:Fi) { # Fi = fish species
      Omega[i,1:C,1:C] ~ dwish(Omega0[1:C,1:C],S)
      pred[i,1:C] ~ dmnorm(mu[i,1:C], Omega[i,,]) #+lenp[i]*lenpred+timep*timepred, Omega[i,,]) # Model prediction.
      for(j in 1:C) { 
        mu[i,j] ~ dnorm(0, 0.0001) # mu1[j], tau1[j]) # Congener-specific mean for fishes
      }
    }
    # Non-dioxins
    for(i in 1:Fi_nd) { # Fi = fish species
      tau_nd[i] ~ dgamma(0.001,0.001)
      for(j in 1:C_nd) {
        pred_[i,j] ~ dnorm(mu[i,j], tau[i])
        mu[i,j] ~ dnorm(0, 0.0001)
      }
    }
  }
  ")

jags <- jags.model(
  mod,
  data = list(
    S = nrow(eu3),
    S_nd = nrow(eu4),
    C = C,
    C_nd = length(conc_nd),
    Fi = Fi,
    Fi_nd = length(fisl_nd),
    conc = log(conc),
    conc_nd = conc_nd,
    #    length = eu3$Length-170, # Subtract average herring size
    #    year = eu3$Year-2009, # Substract baseline year
    fis = match(eu3$Fish, fisl),
    fis_nd = match(eu4$Fish, fisl_nd),
    #    lenpred = 233-170,
    #    timepred = 2009-2009,
    Omega0 = diag(C)/100000
  ),
  n.chains = 4,
  n.adapt = 100
)

update(jags, 100)

samps.j <- jags.samples(
  jags, 
  c(
    'mu', # mean by fish and compound
    'Omega', # precision matrix by compound
    #    'lenp',# parameters for length
    #    'timep', # parameter for Year
    'pred' # predicted concentration for year 2009 and length 17 cm
  ), 
  #  thin=1000,
  N
)
dimnames(samps.j$Omega) <- list(Fish = fisl, Compound = conl, Compound2 = conl, Iter=1:N, Chain=1:4)
dimnames(samps.j$mu) <- list(Fish = fisl, Compound = conl, Iter = 1:N, Chain = 1:4)
#dimnames(samps.j$lenp) <- list(Fish = fisl, Iter = 1:N, Chain = 1:4)
dimnames(samps.j$pred) <- list(Fish = fisl, Compound = conl, Iter = 1:N, Chain = 1:4)
#dimnames(samps.j$timep) <- list(Dummy = "time", Iter = 1:N, Chain = 1:4)

##### conc.param contains expected values of the distribution parameters from the model

conc.param <- Ovariable(
  "conc.param",
  dependencies = data.frame(Name = "samps.j", Ident=NA),
  formula = function(...) {
    conc.param <- list(
      Omega = apply(samps.j$Omega, MARGIN = 1:3, FUN = mean),
      #      lenp = cbind(
      #        mean = apply(samps.j$lenp, MARGIN = 1, FUN = mean),
      #        sd = apply(samps.j$lenp, MARGIN = 1, FUN = sd)
      #      ),
      mu = apply(samps.j$mu, MARGIN = 1:2, FUN = mean)#,
      #      timep = cbind(
      #        mean = apply(samps.j$timep, MARGIN = 1, FUN = mean),
      #        sd = apply(samps.j$timep, MARGIN = 1, FUN = sd)
      #      )
    )
    #    names(dimnames(conc.param$lenp)) <- c("Fish","Metaparam")
    #    names(dimnames(conc.param$timep)) <- c("Dummy","Metaparam")
    conc.param <- melt(conc.param)
    colnames(conc.param)[colnames(conc.param)=="value"] <- "Result"
    colnames(conc.param)[colnames(conc.param)=="L1"] <- "Parameter"
    conc.param$Dummy <- NULL
    #    conc.param$Metaparam <- ifelse(is.na(conc.param$Metaparam), conc.param$Parameter, as.character(conc.param$Metaparam))
    return(Ovariable(output=conc.param, marginal=colnames(conc.param)!="Result"))
  }
)

objects.store(conc.param, samps.j)
cat("Lists conc.params and samps.j stored.\n")

######################3

cat("Descriptive statistics:\n")

# Leave only the main fish species and congeners and remove others

#oprint(summary(
#  eu2[eu2$Compound %in% indices$Compound.PCDDF14 & eu$Fish %in% fisl , ],
#  marginals = c("Fish", "Compound"), # Matrix is always 'Muscle'
#  function_names = c("mean", "sd")
#))

tmp <- eu2[eu2$Compound %in% c("PCDDF","PCB","BDE153","PBB153","PFOA","PFOS","DBT","MBT","TBT"),]@output
ggplot(tmp, aes(x = eu2Result, colour=Fish))+stat_ecdf()+
  facet_wrap( ~ Compound, scales="free_x")+scale_x_log10()

conc.param <- EvalOutput(conc.param)

if(FALSE) {
  scatterplotMatrix(t(exp(samps.j$pred[1,,,1])), main = "Predictions for all compounds for Baltic herring")
  scatterplotMatrix(t(exp(samps.j$pred[,1,,1])), main = "Predictions for all fish species for PCDDF")
  scatterplotMatrix(t(samps.j$Omega[,1,1,,1]))
  #scatterplotMatrix(t(cbind(samps.j$Omega[1,1,1,,1],samps.j$mu[1,1,,1])))
  
  plot(coda.samples(jags, 'Omega', N))
  plot(coda.samples(jags, 'mu', N))
  plot(coda.samples(jags, 'lenp', N))
  plot(coda.samples(jags, 'timep', N))
  plot(coda.samples(jags, 'pred', N))
}