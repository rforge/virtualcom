.transf <- function(x){
	if(min(x)<=0) x <- x - min(x) 
	return (x + 0.0000001)
}
