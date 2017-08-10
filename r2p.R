r2p<-function(r, n) {
	t<-r*sqrt((n-2)/(1-r^2))
	p.value <- 2*pt(-abs(t),n-2)
}