checkBmwpFam <- function(df, famNames, stNames){
  famDf <- df[["Family"]]
  famDf$Family <- as.character(famDf$Family)
  famCheck <- famDf[which(famDf[,"Family"] %in% famNames[,"Family"]),"Family"]
  if(length(famCheck)>0){
    for(i in 1:length(famCheck)){
      taxName <- famCheck[i]
      famDf$Family[famDf$Family == taxName] <- as.character(subset(famNames, Family==taxName)[,2])
    }
    famDf <- aggregate(. ~ Family, famDf, sum)
    names(famDf) <- c("Family", stNames)
    df[["Family"]] <- famDf
    return(df)
  }
  else{df <- df
  return(df)
  }
}
