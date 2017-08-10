fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
pdiffcor<-function(r1, r2, n1, n2){
  statistic <- (fisher.r2z(r1) - fisher.r2z(r2))/sqrt(1/(n1 -3) + 1/(n2 - 3))
  p.value <- 2*pnorm(-abs(statistic))
}
