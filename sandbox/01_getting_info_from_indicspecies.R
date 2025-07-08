# testing ability to get required info from indicspecies

# SETUP ####

# dependencies
library(indicspecies)

# example data
# A data frame with 41 sites (rows) and 33 species (columns).
# Abundance values are represented in abundance classes.
data(wetland)
# define arbitrary groups in a vector
groups <- c(rep(1, 17), rep(2, 14), rep(3,10))




# TEST ANALYSIS ####

# indicator species analysis
# we can think about how to allow implementation of different control and function options (...)
indval <- multipatt(wetland, groups,
                    control = how(nperm=999))




# EXPLORING OUTPUT ####

# investigating the output from multipatt
names(indval)

# A matrix of specificity values (sometimes called positive predictive value) for each species
indval$A
# A matrix of fidelity values (sometimes called sensitivity) for each species
indval$B
# A data frame with test statistics, p-values, and best group combinations for each species
indval$sign
# A list showing, for each species, which combination of groups it is associated with
indval$str
# A matrix showing the group combinations tested. Rows = combinations, columns = groups.
indval$comb

# looks like we want '$sign' ???
out <- indval[['sign']]
out$stat


# MULTIPATT FUNCTION ####
?multipatt
function (x, cluster, func = "IndVal.g", duleg = FALSE, restcomb = NULL,
          min.order = 1, max.order = NULL, control = how(), permutations = NULL,
          print.perm = FALSE)
{
  cl.comb <- function(clnames, min.order, max.order) {
    k <- length(clnames)
    ep <- NULL
    names.ep <- NULL
    for (j in max(min.order, 1):min(max.order, k)) {
      nco <- choose(k, j)
      co <- combn(k, j)
      epn <- matrix(0, ncol = nco, nrow = k)
      for (i in 1:ncol(co)) {
        epn[co[, i], i] <- 1
        if (is.null(names.ep))
          names.ep <- paste(clnames[co[, i]], collapse = "+")
        else names.ep <- c(names.ep, paste(clnames[co[,
                                                      i]], collapse = "+"))
      }
      if (is.null(ep))
        ep <- epn
      else ep <- cbind(ep, epn)
    }
    colnames(ep) <- names.ep
    return(ep)
  }
  rcomb <- function(x, memb, comb, min.order, max.order, mode = "group",
                    restcomb = NULL) {
    k = ncol(memb)
    nsps = ncol(x)
    N = dim(comb)[1]
    ni = colSums(memb * memb)
    tx <- t(x)
    aisp = (tx %*% memb)
    lisp = (tx^2 %*% memb)
    if (mode == "site") {
      lspK = rowSums(lisp)
      aspK = rowSums(aisp)
      aspC = (tx %*% comb)
      nC = colSums(comb * comb)
    }
    else if (mode == "group") {
      aispni = sweep(aisp, 2, ni, "/")
      lispni = sweep(lisp, 2, ni, "/")
      lspK = (N/k) * rowSums(lispni)
      aspK = (N/k) * rowSums(aispni)
      if (max.order == 1) {
        aspC = (N/k) * aispni
        nC = rep(N/k, k)
      }
      else {
        aspC = matrix(0, nrow = nsps, ncol = ncol(comb))
        nC = vector(mode = "numeric", length = ncol(comb))
        cnt = 1
        for (level in min.order:max.order) {
          co = combn(1:k, level)
          for (j in 1:ncol(co)) {
            if (nrow(co) > 1)
              aspC[, cnt] = rowSums(aispni[, co[, j]])
            else aspC[, cnt] = aispni[, co[, j]]
            nC[cnt] = length(co[, j])
            cnt = cnt + 1
          }
        }
        aspC = (N/k) * aspC
        nC = (N/k) * nC
      }
    }
    num = N * aspC - aspK %o% nC
    den = sqrt(((N * lspK) - aspK^2) %o% (N * nC - (nC^2)))
    str = num/den
    colnames(str) <- colnames(comb)[1:ncol(str)]
    if (!is.null(restcomb)) {
      if (sum(restcomb %in% (1:ncol(str))) != length(restcomb))
        stop(paste("One or more indices in 'restcomb' are out of range [1, ",
                   ncol(str), "]", sep = ""))
      str <- str[, restcomb, drop = FALSE]
    }
    return(str)
  }
  indvalcomb <- function(x, memb, comb, min.order, max.order,
                         mode = "group", restcomb = NULL, indvalcomp = FALSE) {
    k <- ncol(memb)
    tx <- t(x)
    aisp = tx %*% comb
    dx <- dim(tx)
    nisp <- matrix(as.logical(tx), nrow = dx[1], ncol = dx[2]) %*%
      comb
    ni = colSums(comb * comb)
    nispni = sweep(nisp, 2, ni, "/")
    if (mode == "site")
      A = sweep(aisp, 1, colSums(x), "/")
    else {
      aispni = sweep(tx %*% memb, 2, colSums(memb), "/")
      asp = rowSums(aispni)
      if (max.order == 1) {
        A <- sweep(aispni, 1, asp, "/")
      }
      else {
        s <- NULL
        for (j in min.order:max.order) {
          if (j == 1) {
            s <- aispni
          }
          else {
            co <- combn(k, j)
            sn <- apply(co, 2, function(x) rowSums(aispni[,
                                                          x]))
            if (is.null(s))
              s <- sn
            else s <- cbind(s, sn)
          }
        }
        A = sweep(s, 1, asp, "/")
      }
    }
    iv = sqrt(A * nispni)
    colnames(iv) <- colnames(comb)
    colnames(A) <- colnames(comb)
    colnames(nispni) <- colnames(comb)
    rownames(A) <- rownames(iv)
    rownames(nispni) <- rownames(iv)
    if (!is.null(restcomb)) {
      if (sum(restcomb %in% (1:ncol(iv))) != length(restcomb))
        stop(paste("One or more indices in 'restcomb' are out of range [1, ",
                   ncol(iv), "]", sep = ""))
      iv = iv[, restcomb, drop = FALSE]
      A = A[, restcomb, drop = FALSE]
      nispni = nispni[, restcomb, drop = FALSE]
    }
    if (!indvalcomp)
      return(iv)
    else return(list(A = A, B = nispni, iv = iv))
  }
  x <- as.matrix(x)
  vegnames <- colnames(x)
  nsps = ncol(x)
  clnames = levels(as.factor(cluster))
  k = length(clnames)
  if (!is.null(control)) {
    nperm = control$nperm
  }
  else if (!is.null(permutations)) {
    nperm = nrow(permutations)
  }
  else {
    stop("You must control permutations with how() or supply a matrix of permutations")
  }
  func = match.arg(func, c("r", "r.g", "IndVal.g", "IndVal",
                           "indval", "indval.g"))
  if (k < 2)
    stop("At least two clusters are required.")
  if (sum(is.na(cluster)) > 0)
    stop("Cannot deal with NA values. Remove and run again.")
  if (sum(is.na(x)) > 0)
    stop("Cannot deal with NA values. Remove and run again.")
  if (!is.null(min.order)) {
    if (!is.numeric(min.order))
      stop("Parameter min.order must be an integer.")
    if (func == "IndVal.g" || func == "indval.g" || func ==
        "IndVal" || func == "indval")
      min.order = min(max(round(min.order), 1), k)
    else min.order = min(max(round(min.order), 1), k - 1)
  }
  else {
    min.order = 1
  }
  if (is.null(max.order)) {
    if (func == "IndVal.g" || func == "indval.g" || func ==
        "IndVal" || func == "indval")
      max.order = k
    else max.order = k - 1
  }
  else {
    if (!is.numeric(max.order))
      stop("Parameter max.order must be an integer.")
    if (func == "IndVal.g" || func == "indval.g" || func ==
        "IndVal" || func == "indval")
      max.order = min(max(round(max.order), 1), k)
    else max.order = min(max(round(max.order), 1), k - 1)
    if (max.order < min.order)
      stop("The value of 'max.order' cannot be smaller than that of 'min.order'.")
  }
  if (!is.null(restcomb)) {
    restcomb = as.integer(restcomb)
  }
  if (duleg)
    max.order = 1
  combin <- cl.comb(clnames, min.order, max.order)
  restcomb = restcomb[restcomb <= ncol(combin)]
  clind = apply(sapply(clnames, "==", cluster), 1, which)
  comb <- combin[clind, ]
  memb <- diag(1, k, k)[clind, ]
  A = NULL
  B = NULL
  if (func == "r")
    str = rcomb(x, memb, comb, min.order, max.order, mode = "site",
                restcomb = restcomb)
  else if (func == "r.g")
    str = rcomb(x, memb, comb, min.order, max.order, mode = "group",
                restcomb = restcomb)
  else if (func == "IndVal" || func == "indval") {
    IndVal = indvalcomb(x, memb, comb, min.order, max.order,
                        mode = "site", restcomb = restcomb, indvalcomp = TRUE)
    str = IndVal$iv
    A = IndVal$A
    B = IndVal$B
  }
  else if (func == "IndVal.g" || func == "indval.g") {
    IndVal = indvalcomb(x, memb, comb, min.order, max.order,
                        mode = "group", restcomb = restcomb, indvalcomp = TRUE)
    str = IndVal$iv
    A = IndVal$A
    B = IndVal$B
  }
  maxstr = apply(str, 1, max)
  wmax <- max.col(str)
  if (!is.null(restcomb))
    m <- as.data.frame(t(combin[, restcomb, drop = FALSE][,
                                                          wmax]))
  else m <- as.data.frame(t(combin[, wmax]))
  row.names(m) <- vegnames
  names(m) <- sapply(clnames, function(x) paste("s", x, sep = "."))
  m$index <- wmax
  m$stat <- apply(str, 1, max)
  pv <- 1
  for (p in 1:nperm) {
    if (!is.null(control)) {
      pInd = shuffle(length(cluster), control = control)
    }
    else {
      pInd = permutations[p, ]
    }
    tmpclind = clind[pInd]
    combp = combin[tmpclind, ]
    membp = diag(1, k, k)[tmpclind, ]
    tmpstr <- switch(func, r = rcomb(x, membp, combp, min.order,
                                     max.order, mode = "site", restcomb = restcomb), r.g = rcomb(x,
                                                                                                 membp, combp, min.order, max.order, mode = "group",
                                                                                                 restcomb = restcomb), indval = indvalcomb(x, membp,
                                                                                                                                           combp, min.order, max.order, mode = "site", restcomb = restcomb),
                     IndVal = indvalcomb(x, membp, combp, min.order, max.order,
                                         mode = "site", restcomb = restcomb), indval.g = indvalcomb(x,
                                                                                                    membp, combp, min.order, max.order, mode = "group",
                                                                                                    restcomb = restcomb), IndVal.g = indvalcomb(x,
                                                                                                                                                membp, combp, min.order, max.order, mode = "group",
                                                                                                                                                restcomb = restcomb))
    tmpmaxstr <- vector(length = nrow(tmpstr))
    for (i in 1:nrow(tmpstr)) tmpmaxstr[i] <- max(tmpstr[i,
    ])
    pv = pv + (tmpmaxstr >= m$stat)
  }
  m$p.value <- pv/(1 + nperm)
  m$p.value[m$index == (2^k - 1)] <- NA
  if (!is.null(restcomb))
    comb <- comb[, restcomb, drop = FALSE]
  a = list(call = match.call(), func = func, cluster = cluster,
           comb = comb, str = str, A = A, B = B, sign = m)
  class(a) = "multipatt"
  return(a)
}
