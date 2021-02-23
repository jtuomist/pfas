# This file contains codes that also exist at Opasnet1836 on page 
# [[Domestic fish consumption of the general population in Finland]]
# The first part was used to preprocess EFSA data and upload it to Opasnet:
# EFSA_fish_consumption_consumers_fi_g_per_day.csv and
# EFSA_fish_consumption_consumers_fi_g_per_day.zip. The zip file did not work so a csv was used.
# They contain the same preprocessed data.


# http://www.efsa.europa.eu/sites/default/files/chronicgdaytotpop.xlsx
efsa1 <- read_xlsx("data/chronicgdaytotpop.xlsx", sheet="L1_All_subjects_g_day", skip=2) # total population
efsa1 <- efsa1[efsa1$Country=="Finland" & grepl("Fish", efsa1$`Foodex L1`),]

efsa3 <- read_xlsx("data/chronicgdaytotpop.xlsx", sheet="L3_All_subjects_g_day", range="R3C1:R131491C16")
efsa3 <- efsa3[efsa3$Country=="Finland" & grepl("Fish", efsa3$`Foodex L2`),]

#tmp <- aggregate(efsa3$Mean, by=efsa3["Foodex L3"], FUN=sum)
#tmp <- tmp[order(-tmp$x),]
#                                            Foodex L3           x
#38                      Salmon and trout (Salmo spp.) 74.26157419
#18                                          Fish meat 62.54528672
#49                                     Tuna (Thunnus) 10.16391067
#50                              Whitefish (Coregonus)  9.87886514
#30                                   Herring (Clupea)  9.22486875
#12                             Fish and potatoes meal  7.58100000
#35                                      Perch (Perca)  4.43950031
#23                                      Fish products  2.38921114
#8                        Cod and whiting (Gadus spp.)  2.00739494
#16                                         Fish balls  1.87321838
#17                                       Fish fingers  1.87172960
#24                                           Fish roe  1.09550296
#26                      Flounder (Platichthys flesus)  0.76270739
#1                                 Anchovy (Engraulis)  0.66095951
#40                        Sardinops (Sardinops sagax)  0.50412277
#42                                Seafood-based meals  0.48000000
#5                                      Bream (Charax)  0.32305245
#33                                 Mackeral (Scomber)  0.26503085
#39                     Sardine and pilchard (Sardina)  0.24428571
#41             Sea catfish and wolf-fish (Anarhichas)  0.15912625
#21                                         Fish paste  0.13883529
#28                                  Hake (Merluccius)  0.06756757
#29                        Halibut (Hippoglossus spp.)  0.03500000
#2                                      Babel (Barbus)  0.00000000
#3                                       Bass (Marone)  0.00000000
# the rest is zero

# http://www.efsa.europa.eu/sites/default/files/chronicgdayconsumers.xlsx
efsa1c <- read_xlsx("data/chronicgdayconsumers.xlsx", sheet="L1_Consumers_only_g_day", range="R3C1:R1132C15") # consumers only
#                    col_types=c(rep("text",5),rep("numeric",10)))
efsa1c <- efsa1c[efsa1c$Country=="Finland" & grepl("Fish", efsa1c$`Foodex L1`),]

efsa3c <- read_xlsx("data/chronicgdayconsumers.xlsx", sheet="L3_Consumers_only_g_day", range="R3C1:R18488C16")
efsa3c <- efsa3c[efsa3c$Country=="Finland" & grepl("Fish", efsa3c$`Foodex L2`),]
efsa3c$nr <- efsa3c$`Nr Consumer` / efsa3c$`% Consumers`
efsa3c$coeff <- efsa3c$STD / efsa3c$Mean
#mean(efsa3c$coeff[efsa3c$coeff>0])
#[1] 0.7353969
#sd(efsa3c$coeff[efsa3c$coeff>0])
#[1] 0.2678654
# It seems that the coefficient of variation is pretty close to 0.75 across all fish species and age groups.
# Therefore, we could make a model where we pool many level3 categories into e.g. imported fish and use local species separately.
# Useful groups: Salmon, fish meat, tuna, herring, whitefish, fish products (incl fish and potatoes meal), 
# perch, flounder, and other species.
# Sample from the ecdf as many as there are actual consumers and form new groups.

out <- data.frame() # This has as many rows as there were consuming individuals on the original studies.
for(i in 1:nrow(efsa3c)) {
  out <- rbind(out, data.frame(
    Row = i,
    Study = efsa3c$Survey[i],
    Participants = efsa3c$nr[i],
    Group = efsa3c$`Pop Class`[i],
    Food = efsa3c$`Foodex L3`[i],
    Result = approx(x=c(0.05, 0.1, 0.5, 0.95, 0.975, 0.99), y=unlist(efsa3c[i,11:16]), n=efsa3c$`Nr Consumer`[i])$y
  ))
}

# write_csv(out, "data/EFSA_fish_consumption_consumers_fi_g_per_day.csv")

tmp <- efsa3c[c(3,5,11:16)] # quantiles saved
colnames(tmp)[colnames(tmp)=="Median"] <- "P50"
tmp <- melt(tmp, id.vars = c("Pop Class","Foodex L3"))
tmp$variable <- as.numeric(gsub("P","",tmp$variable))

ggplot(tmp, aes(x=value, y=variable, colour=`Pop Class`, group=`Pop Class`))+geom_line()+ # construct from data and quantiles
  facet_wrap(~`Foodex L3`, scales="free_x")

ggplot(out, aes(x=Result, color=Group))+stat_ecdf()+ # construct from sampled data
  facet_wrap(~Food, scales="free_x")

# The two graphs above illustrate that the consumer-only quantiles can successfully be used to approximate values and
# create a dataset for further sampling and analysis. Now food groups can be combined before sampling.
# Individuals' consumption is then created by sampling each approximation with the probability given (% Consumers)
# and then summing up to useful categories.

#####################################################3

amount <- Ovariable(
  "amount",
  dependencies = data.frame(Name="dummy"),
  formula = function(...) {
    
    dat <- opasnet.csv("8/88/EFSA_fish_consumption_consumers_fi_g_per_day.csv",
                       wiki = "opasnet_en", sep=",", dec=".", header=TRUE)
    
    out <- data.frame()
    grouping <- unique(dat[c("Group","Food")])
    if(openv$N>0) {
      for(i in 1:nrow(grouping)) {
        tmp <- dat[dat$Group == grouping$Group[i] & dat$Food == grouping$Food[i] , ]
        out <- rbind(out, data.frame(
          Iter = 1:openv$N,
          Group = tmp$Group[1],
          Food = tmp$Food[1],
          Result = sample(tmp$Result, size=openv$N, replace=TRUE) * rbernoulli(openv$N, p=(nrow(tmp) / sum(unique(tmp$Participants)))) # There may be several studies but they never have the exact same number of patients.
        ))
      }
    } else {
      for(i in 1:nrow(grouping)) {
        tmp <- dat[dat$Group == grouping$Group[i] & dat$Food == grouping$Food[i] , ]
        out <- rbind(out, data.frame(
          Group = tmp$Group[1],
          Food = tmp$Food[1],
          Result = median(tmp$Result) * (nrow(tmp) / sum(unique(tmp$Participants)))
        ))
      }
    }
    
    # sort(unique(efsa3c$`Foodex L3`))
    fishes <- data.frame(
      Food = c(
        "Anchovy (Engraulis)",                    "Bream (Charax)",                        
        "Cod and whiting (Gadus spp.)",           "Fish and potatoes meal",                
        "Fish balls",                             "Fish fingers",                          
        "Fish meat",                              "Fish paste",                            
        "Fish products",                          "Fish roe",                              
        "Flounder (Platichthys flesus)",          "Hake (Merluccius)",                     
        "Halibut (Hippoglossus spp.)",            "Herring (Clupea)",                      
        "Mackeral (Scomber)",                     "Perch (Perca)",                         
        "Salmon and trout (Salmo spp.)",          "Sardine and pilchard (Sardina)",        
        "Sardinops (Sardinops sagax)",            "Sea catfish and wolf-fish (Anarhichas)",
        "Seafood-based meals",                    "Tuna (Thunnus)",                        
        "Whitefish (Coregonus)"),
      Fish = c("Imported fish","Bream","Imported fish",rep("Fish product",7),"Flounder",rep("Imported fish",2),
               "Herring","Imported fish", "Perch","Salmon",rep("Imported fish",4), "Tuna","Whitefish")
    )
    
    out <- merge(out, fishes, all.x=TRUE)
    out <- Ovariable(output=out, marginal = colnames(out)!="Result")
    return(out)
  },
  unit="g/d"
)

##############################

samp <- EvalOutput(amount)@output

#samp <- aggregate(samp["amountResult"], by = samp[c("Group","Fish","Iter"[openv$N>0])], FUN=sum)

# mean(samp$amountResult > 0 & samp$amountResult<1)
# [1] 0.000586
# We can cut the x axis from 1 g/day without losing data.

ggplot(samp, aes(x=amountResult+0.1, color=Group))+stat_ecdf()+ # construct from sampled data
  facet_wrap(~Food, scales="free_x")+scale_x_log10()

ggplot(samp, aes(x=amountResult+0.1, color=Group))+stat_ecdf()+ # construct from sampled data
  facet_wrap(~Fish)+scale_x_log10()+coord_cartesian(xlim=c(1,300))+
  labs(
    title="Kalansyönti Suomessa iän ja lajin mukaan",
    y="Kumulatiivinen todennäköisyys",
    x="Kalansyönti (g/päivä)"
  )

tmp <- aggregate(samp["amountResult"], by = samp[c("Group","Iter"[openv$N>0])], FUN=sum)
ggplot(tmp, aes(x=amountResult+0.1, color=Group))+stat_ecdf()+ # construct from sampled data
  scale_x_log10()+coord_cartesian(xlim=c(1,300))+
  labs(
    title="Kokonaiskalansyönti eri ikäryhmissä",
    y="Kumulatiivinen todennäköisyysjakauma",
    x="Kalansyönti (g/päivä"
  )
