zero_dist_traits <- function(x, mat_dissim, BIN) {

  # avoid RCMD notes

  Taxon <- value <- NULL

  trait.df <- as.matrix(mat_dissim)
  diag(trait.df) <- NA
  trait.df[lower.tri(trait.df)] <- NA
  trait.df <- data.frame(Taxon = rownames(trait.df), trait.df)
  trait.df.long <- trait.df %>%
    pivot_longer(-Taxon) %>%
    filter(value == 0) %>%
    mutate(Taxon = as.character(Taxon))

  for (i in 1:nrow(trait.df.long)) {
    name1 <- trait.df.long[i, 1]
    name2 <- trait.df.long[i, 2]
    x[x[, "Taxon"] %in% name2, "Taxon"] <- name1
    trait.df <- trait.df[, !colnames(trait.df) %in% name2]
    trait.df <- trait.df[!rownames(trait.df) %in% name2, ]
  }

  x <- aggregate(. ~ Taxon, x, FUN = sum)

  if (BIN) {
    x <- to_bin(x)
  }

  trait.df <- as.matrix(trait.df[, -1])

  if (any(x[, "Taxon"] != rownames(trait.df))) {
    stop("Something went wrong with the zero distance removal. Please ask the mainteiner to solve the problem.")
  }
  diag(trait.df) <- rep(0, ncol(trait.df))
  trait.df[lower.tri(trait.df)] <- trait.df[upper.tri(trait.df)]
  colnames(trait.df) <- rownames(trait.df)
  trait.dist <- as.dist(trait.df)

  # returns modified Df and distances togheter with the list of taxa with the same traits

  list(x, trait.dist, as.data.frame(trait.df.long[, 1:2]))
}
