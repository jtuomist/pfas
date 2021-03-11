#This is code Op_en3104/conc_poll on page [[EU-kalat]]

library(OpasnetUtils)

#objects.latest("Op_en3104", code_name="pollutant_bayes")

conc_poll <- Ovariable(
  "conc_poll",
  dependencies = data.frame(
    Name=c("conc_param"), #,"lengt","time"),
    Ident=c("Op_en3104/pollutant_bayes")#,NA,NA)
  ),
  formula=function(...) {
    require(MASS)
    tmp1 <- conc_param + Ovariable(data=data.frame(Result="0-1")) # Ensures Iter #  lengt + time + 
    tmp2 <- unique(tmp1@output[setdiff(
      colnames(tmp1@output)[tmp1@marginal],
      c("Compound","Compound2","Metaparam","Parameter")
    )])
    tmp2$Row <- 1:nrow(tmp2)
    tmp3 <- merge(tmp2,tmp1@output)
    out <- data.frame()
    for(i in 1:nrow(tmp2)) {
      tmp <- tmp3[tmp3$Row == i , ]
      Omega <- solve(tapply( # Is it sure that PCDDF and PCB are not mixed to wrong order?
        tmp$conc_paramResult[tmp$Parameter=="Omega"],
        tmp[tmp$Parameter=="Omega", c("Compound","Compound2")],
        sum # Equal to identity because only 1 row per cell.
      )) # Precision matrix
      con <- names(Omega[,1])
      
      mu <- tmp$conc_paramResult[tmp$Parameter=="mu"][match(con,tmp$Compound[tmp$Parameter=="mu"])] # + # baseline
#        rnorm(1,
#              tmp$conc_paramResult[tmp$Parameter=="lenp" & tmp$Metaparam=="mean"][1],
#              tmp$conc_paramResult[tmp$Parameter=="lenp" & tmp$Metaparam=="sd"][1]
#        ) * (tmp$lengtResult[1]-170) + # lengt
#        rnorm(1,
#              tmp$conc_paramResult[tmp$Parameter=="timep" & tmp$Metaparam=="mean"][1],
#              tmp$conc_paramResult[tmp$Parameter=="timep" & tmp$Metaparam=="sd"][1]
#        )* (tmp$timeResult[1]-2009) # time
      
      rnd <- exp(mvrnorm(1, mu, Omega))
      out <- rbind(out, merge(tmp2[i,], data.frame(Compound=con,Result=rnd)))
    }
    out$Row <- NULL
#    temp <- aggregate(
#      out["Result"],
#      by=out[setdiff(colnames(out), c("Result","Compound"))],
#      FUN=sum
#    )
#    temp$Compound <- "TEQ"
    out <- Ovariable(
      output = out, # rbind(out, temp),
      marginal = colnames(out) != "Result"
    )
    return(out)
  }
)

objects.store(conc_poll)
cat("Ovariable conc_poll stored.\n")
