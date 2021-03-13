# This is code Op_en3104/pollutant_bayes on page [[EU-kalat]]
# The code is also available at https://github.com/jtuomist/pfas/blob/main/conc_pcddf_preprocess.R

library(OpasnetUtils)
library(reshape2)
library(rjags) # JAGS
library(ggplot2)
library(MASS) # mvrnorm
library(car) # scatterplotMatrix

#' Find the level of quantification for dinterval function
#' @param df data.frame
#' @return data.matrix
add_loq <- function(df) { # This should reflect the fraction of observations below LOQ.
  LOQ <- unlist(lapply(df, FUN = function(x) min(x[x!=0], na.rm=TRUE))) 
  out <- sapply(
    1:length(LOQ),
    FUN = function(x) ifelse(df[,x]==0, 0.5*LOQ[x], df[,x])
  )
  out <- data.matrix(out)
  return(out)
}

#size <- Ovariable("size", ddata="Op_en7748", subset="Size distribution of fish species")
#time <- Ovariable("time", data = data.frame(Result=2015))

objects.latest("Op_en3104", code_name = "preprocess2") # [[EU-kalat]] euw

# Hierarchical Bayes model.

# PCDD/F concentrations in fish.
# It uses the TEQ sum of PCDD/F (PCDDF) as the total concentration
# of dioxin and PCB respectively for PCB in fish.
# PCDDF depends on size of fish, fish species, catchment time, and catchment area,
# but we omit catchment area. In addition, we assume that size of fish has
# zero effect for other fish than Baltic herring.
# Catchment year affects all species similarly. 

eu3 <- euw[!colnames(euw) %in% c("MPhT","DOT","BDE138")] # No values > 0

eu3 <- eu3[eu3$Matrix == "Muscle" , ]
eu3$Locat <- ifelse(eu3$Location=="Porvoo",2,
                       ifelse(eu3$Location=="Helsinki, Vanhankaupunginlahti Bay",3,1))
locl <- c("Finland","Porvoo","Helsinki")

#conl_nd <- c("PFAS","PFOA","PFOS","DBT","MBT","TBT","DPhT","TPhT")
conl_nd <- c("PFAS","PFOS") # TBT would drop Porvoo measurements
fisl <- fisl_nd <- c("Baltic herring","Bream","Flounder","Perch","Roach","Salmon","Whitefish")

eu4 <- eu3[rowSums(is.na(eu3[conl_nd]))<length(conl_nd) & eu3$Fish %in% fisl_nd ,
           c(1:5,match(c("Locat",conl_nd),colnames(eu3)))]

conc_nd <- add_loq(eu4[conl_nd])

conl <- c("TEQ","PCDDF","PCB") # setdiff(colnames(eu3)[-(1:5)], conl_nd)
eu3 <- eu3[!is.na(eu3$PCDDF) & eu3$Fish %in% fisl , c(1:5, match(conl,colnames(eu3)))]

oprint(head(eu3))
oprint(head(eu4))

C <- length(conl)
Fi <- length(fisl)
N <- 200
thin <- 100
conl
fisl
conl_nd
fisl_nd

eu3 <- eu3[rowSums(is.na(eu3))==0,]
conc <- add_loq(eu3[conl]) # Remove rows with missing data.

# The model assumes that all fish groups have the same Omega but mu varies.

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
        conc_nd[i,j] ~ dnorm(muind_nd[i,j], tau_nd[j])
        muind_nd[i,j] <- mu_nd[fis_nd[i],j] + mulocat[locat[i]] #+ lenp[fis[i]]*length[i] + timep*year[i]
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
    mulocat[1] <- 0
    mulocat[2] ~ dnorm(0,0.001)
    mulocat[3] ~ dnorm(0,0.001)
    for(j in 1:C_nd) {
      tau_nd[j] ~ dgamma(0.001,0.001)
      for(i in 1:Fi_nd) { # Fi = fish species
        pred_nd[i,j] ~ dnorm(mu[i,j], tau_nd[j])
        mu_nd[i,j] ~ dnorm(0, 0.0001)
      }
    }
  }
")

jags <- jags.model(
  mod,
  data = list(
    S = nrow(conc),
    S_nd = nrow(conc_nd),
    C = C,
    C_nd = ncol(conc_nd),
    Fi = Fi,
    Fi_nd = length(fisl_nd),
    conc = log(conc),
    conc_nd = log(conc_nd),
    locat = eu4$Locat,
    #    length = eu3$Length-170, # Subtract average herring size
    #    year = eu3$Year-2009, # Substract baseline year
    fis = match(eu3$Fish, fisl),
    fis_nd = match(eu4$Fish, fisl_nd),
    #    lenpred = 233-170,
    #    timepred = 2009-2009,
    Omega0 = diag(C)/100000
  ),
  n.chains = 4,
  n.adapt = 200
)

update(jags, 1000)

samps.j <- jags.samples(
  jags, 
  c(
    'mu', # mean by fish and compound
    'Omega', # precision matrix by compound
    #    'lenp',# parameters for length
    #    'timep', # parameter for Year
    'pred', # predicted concentration for year 2009 and length 17 cm
    'pred_nd',
    'mu_nd',
    'tau_nd',
    'mulocat'
  ), 
  thin=thin,
  N*thin
)
dimnames(samps.j$Omega) <- list(Fish = fisl, Compound = conl, Compound2 = conl, Iter=1:N, Chain=1:4)
dimnames(samps.j$mu) <- list(Fish = fisl, Compound = conl, Iter = 1:N, Chain = 1:4)
#dimnames(samps.j$lenp) <- list(Fish = fisl, Iter = 1:N, Chain = 1:4)
dimnames(samps.j$pred) <- list(Fish = fisl, Compound = conl, Iter = 1:N, Chain = 1:4)
dimnames(samps.j$pred_nd) <- list(Fish = fisl_nd, Compound = conl_nd, Iter = 1:N, Chain = 1:4)
dimnames(samps.j$mu_nd) <- list(Fish = fisl_nd, Compound = conl_nd, Iter = 1:N, Chain = 1:4)
dimnames(samps.j$tau_nd) <- list(Compound = conl_nd, Iter = 1:N, Chain = 1:4)
#dimnames(samps.j$timep) <- list(Dummy = "time", Iter = 1:N, Chain = 1:4)
dimnames(samps.j$mulocat) <- list(Area = locl, Iter = 1:N, Chain = 1:4)

##### conc_param contains expected values of the distribution parameters from the model

conc_param <- list(
  Omega = apply(samps.j$Omega, MARGIN = 1:3, FUN = mean),
  #      lenp = cbind(
  #        mean = apply(samps.j$lenp, MARGIN = 1, FUN = mean),
  #        sd = apply(samps.j$lenp, MARGIN = 1, FUN = sd)
  #      ),
  mu = apply(samps.j$mu, MARGIN = 1:2, FUN = mean),
  #      timep = cbind(
  #        mean = apply(samps.j$timep, MARGIN = 1, FUN = mean),
  #        sd = apply(samps.j$timep, MARGIN = 1, FUN = sd)
  #      )
  mu_nd =  apply(samps.j$mu_nd, MARGIN = 1:2, FUN = mean),
  tau_nd =  apply(samps.j$tau_nd, MARGIN = 1, FUN = mean),
  mulocat = apply(samps.j$mulocat, MARGIN = 1, FUN = mean)
)
#    names(dimnames(conc_param$lenp)) <- c("Fish","Metaparam")
#    names(dimnames(conc_param$timep)) <- c("Dummy","Metaparam")

conc_param <- melt(conc_param)
colnames(conc_param)[colnames(conc_param)=="value"] <- "Result"
colnames(conc_param)[colnames(conc_param)=="L1"] <- "Parameter"
conc_param$Compound[conc_param$Parameter =="tau_nd"] <- conl_nd # drops out because one-dimensional
conc_param$Area[conc_param$Parameter =="mulocat"] <- locl # drops out because one-dimensional
conc_param <- fillna(conc_param,c("Fish","Area"))
for(i in 1:ncol(conc_param)) {
  if("factor" %in% class(conc_param[[i]])) conc_param[[i]] <- as.character(conc_param[[i]])
}
conc_param <- Ovariable("conc_param",data=conc_param)

objects.store(conc_param)
cat("Ovariable conc_param stored.\n")

######################3

cat("Descriptive statistics:\n")

# Leave only the main fish species and congeners and remove others

#oprint(summary(
#  eu2[eu2$Compound %in% indices$Compound.PCDDF14 & eu$Fish %in% fisl , ],
#  marginals = c("Fish", "Compound"), # Matrix is always 'Muscle'
#  function_names = c("mean", "sd")
#))

#tmp <- euw[euw$Compound %in% c("PCDDF","PCB","BDE153","PBB153","PFOA","PFOS","DBT","MBT","TBT"),]
#ggplot(tmp, aes(x = eu2Result, colour=Fish))+stat_ecdf()+
#  facet_wrap( ~ Compound, scales="free_x")+scale_x_log10()

scatterplotMatrix(t(exp(samps.j$pred[2,,,1])), main = paste("Predictions for several compounds for",
                                                            names(samps.j$pred[,1,1,1])[2]))
scatterplotMatrix(t(exp(samps.j$pred[,1,,1])), main = paste("Predictions for all fish species for",
                                                            names(samps.j$pred[1,,1,1])[1]))
scatterplotMatrix(t(samps.j$Omega[2,,1,,1]), main = "Omega for several compounds in Baltic herring")

scatterplotMatrix(t((samps.j$pred_nd[1,,,1])), main = paste("Predictions for several compounds for",
                                                            names(samps.j$pred_nd[,1,1,1])[1]))

scatterplotMatrix(t((samps.j$mulocat[,,1])), main = paste("Predictions for location average difference",
                                                            names(samps.j$pred_nd[,1,1,1])[1]))

#plot(coda.samples(jags, 'Omega', N))
plot(coda.samples(jags, 'mu', N*thin, thin))
#plot(coda.samples(jags, 'lenp', N))
#plot(coda.samples(jags, 'timep', N))
plot(coda.samples(jags, 'pred', N*thin, thin))
plot(coda.samples(jags, 'mu_nd', N*thin, thin))
plot(coda.samples(jags, 'mulocat', N*thin, thin))
tst <- (coda.samples(jags, 'pred', N))
