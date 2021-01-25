manage_inconsistencies <- function(DF, Tree) {
  taxa <- as.character(DF[, "Taxon"])
  taxa <- trimws(taxa)
  taxa <- unlist(lapply(as.list(taxa), capWords))

  # Position of taxon in the df data.frame
  taxind <- data.frame(row = numeric(), col = numeric())
  for (i in 1:length(taxa)) {
    temp <- which(Tree == taxa[i], arr.ind = T)
    taxind <- rbind(temp, taxind)
  }

  # find the Taxon corresponding to taxind and add it to taxind
  taxa.t <- c()
  for (i in 1:nrow(taxind)) {
    taxa.t <- c(taxa.t, as.character(Tree[taxind[i, 1], taxind[i, 2]]))
  }
  taxind <- cbind(taxind, taxa.t)

  # order taxind by column to get finest taxonomic levels up
  taxind <- taxind[order(taxind[, 2], decreasing = FALSE), ]

  # if a row is duplicated there is the need to solve inconsistencies. this is due to the
  # characteristics of the Tree in biomonitoR.

  # count the number of times a row appears and select only those appearing more than once
  r.count <- table(taxind[, 1])
  r.count <- as.numeric(names(r.count[r.count > 1]))
  inconsistency <- taxind[taxind[, 1] %in% r.count, ]

  if (nrow(inconsistency) == 0) {
    return(DF)
  }

  # we work on a copy in order to keep the orginal abundances for further calculations
  DF.temp <- DF

  mes <- "Parent-child pairs found, please see traceB for further information"

  store.taxa <- data.frame(tax_1 = numeric(), tax2 = numeric())
  tax.store <- c()
  for (i in 1:nrow(inconsistency)) {
    temp <- inconsistency[inconsistency[, 1] == inconsistency[i, 1], ]
    temp <- temp[(nrow(temp) - 1):nrow(temp), ]
    Col1 <- temp[1, 2]
    Col2 <- temp[2, 2]
    Row <- inconsistency[i, 1]
    tax2 <- as.character(Tree[Row, Col1])
    tax1 <- as.character(Tree[Row, Col2])
    if (tax1 %in% tax.store & tax2 %in% tax.store) {
      next()
    }
    to.rem1 <- DF.temp[DF.temp$Taxon == tax2, -1, drop = FALSE]
    to.rem2 <- DF.temp[DF.temp$Taxon == tax1, -1, drop = FALSE]
    DF.temp[DF.temp$Taxon == tax2, -1] <- to.rem1 - to.rem2
    store.taxa <- rbind(store.taxa, data.frame(tax_1 = tax1, tax_2 = tax2))
    tax.store <- unique(c(tax.store, tax1, tax2))
  }

  colnames(store.taxa) <- c("Child", "Parent")
  message(mes)
  list(DF = DF.temp, Changes = store.taxa)
}


manage_exceptions <- function(DF, Tree, y, Taxon) {
  taxa <- as.character(Taxon)
  taxa <- trimws(taxa)
  taxa <- unlist(lapply(as.list(taxa), capWords))

  # Position of taxon in the df data.frame
  taxind <- data.frame(row = numeric(), col = numeric())
  for (i in 1:length(taxa)) {
    temp <- which(Tree == taxa[i], arr.ind = T)
    taxind <- rbind(temp, taxind)
  }

  if (nrow(taxind) == 0) {
    return(DF)
  }

  # find the Taxon corresponding to taxind and add it to taxind
  taxa.t <- c()
  for (i in 1:nrow(taxind)) {
    taxa.t <- c(taxa.t, as.character(Tree[taxind[i, 1], taxind[i, 2]]))
  }

  taxind$taxa.t <- as.character(taxa.t)
  taxind <- aggregate(. ~ col + taxa.t, taxind, FUN = min)
  taxind <- taxind[, c(3, 1, 2)]


  # order taxind by column to get finest taxonomic levels up
  taxind <- taxind[order(taxind[, 2], decreasing = FALSE), ]
  taxind <- taxind[!duplicated(taxind[, "taxa.t", drop = FALSE]), ]



  higher.taxa <- list()
  lower.taxa <- list()
  for (i in 1:nrow(taxind)) {
    higher.t <- Tree[taxind[i, 1], 1:(taxind[i, 2])]
    higher.t <- t(higher.t)[, 1]
    higher.taxa[[i]] <- higher.t[higher.t %in% y[, "Taxon"]]
    lower.t <- Tree[Tree[, taxind[i, 2]] == taxind[i, 3], ]
    lower.t <- lower.t[, (taxind[i, 2] + 1):10]
    lower.t <- as.character(unlist(lower.t))
    lower.taxa[[i]] <- lower.t[lower.t %in% y[, "Taxon"]]
  }

  DF.temp <- DF

  mes <- "Exceptions found, please see traceB for further information"

  for (i in 1:length(higher.taxa)) {
    tax <- taxind[i, 3]
    tax.up <- higher.taxa[[i]]
    tax.lo <- lower.taxa[[i]]
    if (length(tax.up) > 0) {
      to.rem <- DF.temp[DF.temp$Taxon %in% tax, -1, drop = FALSE]
      temp.up <- sweep(DF.temp[DF.temp$Taxon %in% tax.up, -1, drop = FALSE], 2, t(to.rem)[, 1], FUN = "-")
      DF.temp[DF.temp$Taxon %in% tax.up, -1] <- temp.up
    }

    if (length(tax.lo) > 0) {
      DF.temp[DF.temp$Taxon %in% tax.lo, -1] <- 0
    }
  }

  store.taxa <- data.frame(tax_1 = c(), tax_2 = c(), operation = c())

  for (i in 1:nrow(taxind)) {
    if (length(higher.taxa[[i]]) > 0) {
      DF1 <- data.frame(tax_1 = rep(taxind[i, "taxa.t"], length(higher.taxa[[i]])), tax_2 = higher.taxa[[i]], operation = rep("removed from", length(higher.taxa[[i]])))
    }
    if (length(lower.taxa[[i]]) > 0) {
      DF2 <- data.frame(tax_1 = rep(taxind[i, "taxa.t"], length(lower.taxa[[i]])), tax_2 = lower.taxa[[i]], operation = rep("set to zero", length(lower.taxa[[i]])))
    }

    if (exists("DF1", inherits = FALSE)) {
      store.taxa <- rbind(store.taxa, DF1)
    }

    if (exists("DF2", inherits = FALSE)) {
      store.taxa <- rbind(store.taxa, DF2)
    }
  }

  rownames(store.taxa) <- NULL
  store.taxa <- store.taxa[!duplicated(store.taxa), ]
  store.taxa <- store.taxa[!store.taxa[, "tax_1"] == store.taxa[, "tax_2"], ]

  message(mes)
  list(DF = DF.temp, Changes = store.taxa[, c("tax_1", "operation", "tax_2")])
}
