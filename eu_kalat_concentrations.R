# This is code Op_fi5932/conc on page [[PFAS-yhdisteiden tautitaakka]]

library(OpasnetUtils)
library(tidyverse)

df <- opbase.data("Op_en3104", subset="POPs") # [[EU-kalat]]
df <- df[df$POP %in% c("PFOA", "PFOS"),-c(8,9,11,14)]
colnames(df) <- c("Code","Matrix","POP","Species","Site","Location","Date","N","Length","Weight","Age","Fat","Dry_matter","Result")
for(i in 8:14) df[[i]] <- as.numeric(df[[i]])
df$Code <- gsub("(M|L)", "", df$Code)
df$Result <- df$Result / 1000 # pg/g f.w --> ng/g f.w.
df$Result[df$Result==0] <- "<0.15" # Level of quantitation

openv.setN(1) # Only one random value as sample
df <- EvalOutput(Ovariable("df",data=df))@output
colnames(df)[colnames(df)=="dfResult"] <- "Result"

ggplot(df, aes(x=Species, y = Result, colour=Location, shape=Matrix))+geom_point()+
  coord_flip()+scale_y_log10()+
  facet_grid(~POP, scales="free_x")+
  labs(
    title="PFAS concentrations in fish in Finland",
    y = "PFAS concentration (ng/g fresh weight)"
  )

aggregate(df$Result, by=df[c("Matrix","Location","Species","POP")], FUN=mean)
# PFOA concentrations are 0 for all except Baltic herring. Therefore, PFOA could
# be removed from further scrutiny.
# Also, remove farmed rainbow trout because it is not representative of its environment.
df <- df[df$POP!="PFOA" & df$Site!="Farmed",]

aggregate(df$Result, by=df[c("Matrix","Location","Species","POP")], FUN=length)

ggplot(df, aes(x=Species, y = Result, colour=Location, shape=Site))+geom_point(size=3)+
  coord_flip()+scale_y_log10()+
  facet_grid(Matrix~., scales="free_x")+
  labs(
    title="PFOS concentrations in Finnish fishes",
    y = "PFOS concentration (ng/g fresh weight)"
  )

v1 <- EvalOutput(Ovariable("v1", data=df[df$Matrix=="Muscle", c(1,3,4,6,14)]))
v2 <- EvalOutput(Ovariable("v2", data=df[df$Matrix=="Liver", c(1,3,4,6,14)]))
tmp <- v2/v1

# It seems that
# 1) Vanhankaupunginlahti in Helsinki seems to be badly polluted, but otherwise minor differences 
#    between places and species in PFOS concentrations. There is some kind of interaction,
#    because pike-perch muscle show similar concentration in Vanhankaupunginlahti and elsewhere,
#    while in perch that place is much more polluted with PFOA, and this is the case with liver as well.
# 2) The liver concentrations seem to be 3 - 30 times higher than in muscle within the same
#    pooled sample. However, there are liver data only about three species: salmon, perch,
#    and burbot. Therefore, liver brings little added value and are not considered later.

tmp <- v1[v1$Location=="Helsinki, Vanhankaupunginlahti Bay" & v1$POP=="PFOS",]
tmp$Iter <- 1:nrow(tmp@output)
summary(tmp, marginals="POP")

v1 <- v1[v1$Location != "Helsinki, Vanhankaupunginlahti Bay",]

ggplot(v1@output, aes(x=v1Result))+stat_ecdf()+
  labs(title="PFAS concentration if fish in Finland (excl Helsinki)",
       x="Concentration (ng/g fresh weight)",
       y="Cumulative probability")

v1$Iter <- 1:length(result(v1)) # Make data probabilistic as if each measurement is a random sample of a single distribution
summary(v1, marginals="POP")
#       POP     mean       sd      min   Q0.025   median   Q0.975      max
# 2.5% PFOS 2055.757 1404.045 305.0399 330.1365 1533.269 5029.697 5814.935

# These concentrations must be in pg/g (fresh weight?)

# Conclusion: All fish species seem to be at the same concentration range, so that is used for all
# fish species. Vanhankaupunginlahti data is excluded because that seems to be polluted like no other place.

# It could be fun to draw PFAS concentrations on map using bubble pie charts. See
# https://stat.ethz.ch/pipermail/r-help/2014-June/375749.html
# https://cran.r-project.org/web/packages/plotrix/plotrix.pdf

#reg <- lm(Result ~ Location + Species, data = df[df$Matrix=="Muscle",])
reg <- lm(Result ~ Species + Site, data = df[df$Matrix=="Muscle",])
summary(reg)

# Based on regression analysis, we can see that there is little difference between
# other sea areas than Helsinki Vanhankaupunginlahti, but clear difference between Sites.
# So, we can drop Location and keep Site. In addition, only two species show up:
# perch have higher PFOS concentrations and pike-perch lower, but otherwise there are
# no statistically significant differences between species.
# So, we could reclassify Species to perch, pike-perch, and other.
