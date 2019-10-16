check_vcov <- function(x){
  if (!isTRUE(all.equal(x, t(x))) || any(diag(x) < 0)){
    saveRDS(x, "c:/git_repositories/BFpack/erordump.RData")
    stop(sQuote("sigma"), " is not a covariance matrix")
  }
}
