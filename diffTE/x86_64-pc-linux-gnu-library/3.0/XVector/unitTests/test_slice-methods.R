test_XDouble_slice <- function() {
  ## Use slice against an Rle as an easy test
  x <- c(0.2, 0.5, 1, 1, 1, 1.5, 1.5, -.5, -.5, -.5, 10.2, 10.3)
  r <- Rle(x)
  
  for (lower in c(-0.5, 0, 1.2, 5)) {
    double.slice <- slice(x, lower)
    rle.slice <- slice(r, lower)
    checkEquals(length(double.slice), length(rle.slice))
    is.same <- sapply(1:length(double.slice), function(i) {
      d <- as.numeric(double.slice[[i]])
      r <- as.numeric(rle.slice[[i]])
      checkEqualsNumeric(d, r)
    })
    checkTrue(all(is.same))
  }
}

