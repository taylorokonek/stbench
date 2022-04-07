# expit
expit <- function(x) {
  exp(x)/(1 + exp(x))
}

# logit
logit <- function(x) {
  log(x/(1-x))
}