# This is code Op_en3104/preprocess2 on page [[EU-kalat]]
library(OpasnetUtils)
library(ggplot2)
library(reshape2)

openv.setN(1)
opts = options(stringsAsFactors = FALSE)
euRaw <- Ovariable("euRaw", ddata = "Op_en3104", subset = "POPs") # [[EU-kalat]]

eu <- Ovariable(
  "eu",
  dependencies = data.frame(
    Name=c("euRaw", "TEF"),
    Ident=c(NA,"Op_en4017/initiate")
  ),
  formula = function(...) {
    out <- euRaw
    out$Length<-as.numeric(as.character(out$Length_mean_mm))
    out$Year <- as.numeric(substr(out$Catch_date, nchar(as.character(out$Catch_date))-3,100))
    out$Weight<-as.numeric(as.character(out$Weight_mean_g))
    
    out <- out[,c(1:6, 8: 10, 14:17, 19:22, 18)] # See below 
    
    #[1] "ﮮTHL_code"             "Matrix"                "POP"                   "Fish_species"         
    #[5] "Catch_site"            "Catch_location"        "Catch_season"          "Catch_square"         
    #[9] "N_individuals"         "Sex"                   "Age"                   "Fat_percentage"       
    #[13] "Dry_matter_percentage" "euRawSource"           "Length"                "Year"                 
    #[17] "Weight"                "euRawResult"          
    
    colnames(out@output)[1:13] <- c("THLcode", "Matrix", "Compound", "Fish", "Site", "Location", "Season",
                                   "Square","N","Sex","Age","Fat","Dry_matter")
    out@marginal <- colnames(out)!="euRawResult"
    
    tmp <- oapply(out * TEF, cols = "Compound", FUN = "sum")
    colnames(tmp@output)[colnames(tmp@output)=="Group"] <- "Compound"
    # levels(tmp$Compound)
    # [1] "Chlorinated dibenzo-p-dioxins" "Chlorinated dibenzofurans"     "Mono-ortho-substituted PCBs"  
    # [4] "Non-ortho-substituted PCBs"   
    levels(tmp$Compound) <- c("PCDD","PCDF","moPCB","noPCB")
    
    out <- OpasnetUtils::combine(out, tmp)
    
    out$Compound <- factor( # Compound levels are ordered based on the data table on [[TEF]]
      out$Compound, 
      levels = unique(c(levels(TEF$Compound), unique(out$Compound)))
    )
    out$Compound <- out$Compound[,drop=TRUE]
    
    return(out)
  }
)

eu <- EvalOutput(eu)

euw <- reshape( 
  eu@output, 
  v.names = "euResult", 
  idvar = c("THLcode", "Matrix", "Fish"), # , "Site","Location","Season","Square","N","Sex","Age","Fat", "Dry_matter","Length","Year","Weight"
  timevar = "Compound", 
  drop = c("euRawSource","TEFversion","TEFrawSource","TEFSource","Source","euSource"), 
  direction = "wide"
)
colnames(euw) <- gsub("euResult\\.","",colnames(euw))
euw$PCDDF <- euw$PCDD + euw$PCDF
euw$PCB <- euw$noPCB + euw$moPCB
euw$TEQ <- euw$PCDDF + euw$PCB
euw$PFAS <- euw$PFOA + euw$PFOS

#################### PFAS measurements from Porvoo

conc_pfas_raw <- EvalOutput(Ovariable(
  "conc_pfas_raw",
  data=opbase.data("Op_fi5932", subset="PFAS concentrations"), # [[PFAS-yhdisteiden tautitaakka]]
  unit="ng/g f.w.")
)@output

conc_pfas_raw <- reshape(conc_pfas_raw,
                         v.names="conc_pfas_rawResult",
                         timevar="Compound",
                         idvar=c("Obs","Fish"),
                         drop="conc_pfas_rawSource",
                         direction="wide")

colnames(conc_pfas_raw) <- gsub("conc_pfas_rawResult\\.","",colnames(conc_pfas_raw))
conc_pfas_raw <- within(conc_pfas_raw, PFAS <- PFOS + PFHxS + PFOA + PFNA)
conc_pfas_raw$Obs <- NULL

euw <- orbind(euw, conc_pfas_raw)

objects.store(euw)
cat("Data.frame euw stored.\n")
