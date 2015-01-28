test_XDoubleViews_equality <- function() {
  x <- rnorm(100)
  bounds <- IRanges(c(1, 20, 50, 80), width=c(5, 10, 15, 18))
  bounds2 <- IRanges(c(10, 30, 50, 80), width=c(5, 8, 15, 18))
  v <- Views(x, bounds)
  v2 <- Views(x, bounds2)
  
  checkTrue(all(v == v))
  checkTrue(all((v != v2) == c(TRUE, TRUE, FALSE, FALSE)))
}


test_XDoubleViews_viewMins <- function() {
  x <- rnorm(100)
  bounds <- IRanges(c(1, 20, 50, 80), width=c(5, 10, 15, 18))
  v <- Views(x, bounds)
  
  val <- viewMins(v)
  expect <- sapply(1:length(bounds), function(i) {
    min(x[start(bounds)[i]:end(bounds[i])])
  })
  
  checkEqualsNumeric(val, expect)
}

test_XDoubleViews_viewMaxs <- function() {
  x <- rnorm(100)
  bounds <- IRanges(c(1, 20, 50, 80), width=c(5, 10, 15, 18))
  v <- Views(x, bounds)
  
  val <- viewMaxs(v)
  expect <- sapply(1:length(bounds), function(i) {
    max(x[start(bounds)[i]:end(bounds[i])])
  })
  checkEqualsNumeric(val, expect)
}

test_XDoubleViews_viewSums <- function() {
  x <- rnorm(100)
  bounds <- IRanges(c(1, 20, 50, 80), width=c(5, 10, 15, 18))
  v <- Views(x, bounds)
  
  val <- viewSums(v)
  expect <- sapply(1:length(bounds), function(i) {
    sum(x[start(bounds)[i]:end(bounds[i])])
  })
  
  checkEqualsNumeric(val, expect)
}

test_XDoubleViews_viewMeans <- function() {
  x <- rnorm(100)
  bounds <- IRanges(c(1, 20, 50, 80), width=c(5, 10, 15, 18))
  v <- Views(x, bounds)
  
  val <- viewMeans(v)
  expect <- sapply(1:length(bounds), function(i) {
    mean(x[start(bounds)[i]:end(bounds[i])])
  })
  
  checkEqualsNumeric(val, expect)
}

test_XDoubleViews_viewWhichMins <- function() {
  x <- rnorm(100)
  bounds <- IRanges(c(1, 20, 50, 80), width=c(5, 10, 15, 18))
  v <- Views(x, bounds)
  
  val <- viewWhichMins(v)
  expect <- sapply(1:length(bounds), function(i) {
    which.min(x[start(bounds)[i]:end(bounds[i])]) + start(bounds)[i] - 1L
  })
  
  checkIdentical(val, expect)
}

test_XDoubleViews_viewWhichMaxs <- function() {
  x <- rnorm(100)
  bounds <- IRanges(c(1, 20, 50, 80), width=c(5, 10, 15, 18))
  v <- Views(x, bounds)
  
  val <- viewWhichMaxs(v)
  expect <- sapply(1:length(bounds), function(i) {
    which.max(x[start(bounds)[i]:end(bounds[i])]) + start(bounds)[i] - 1L
  })
  
  checkIdentical(val, expect)
}
