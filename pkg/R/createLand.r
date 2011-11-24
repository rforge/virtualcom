#-----------------------------------------------------------------------------
# Create spatially correlated landscapes with random fields
#-------------------------------------------------------------------------------
#gridSize:  x and y range of the grid (ex. 0:(10 - 1))
# minEnv=min env value; maxEnv=max env value, meanEnv= mean env value, varEnv=var of env value)
#ltype: random landscape:"lr",  spatially auto-correlated landscape:"la"; gradient landscape:"lg"

library("RandomFields")

createLand <- function(ltype, gridSize, minEnv=0, maxEnv=100, meanEnv=50, varEnv=200 ){
   gridside<- max(gridSize)+1
   if (ltype == "la"){
      ###### spatially auto-correlated landscape
      model <- "stable"
      mean <- meanEnv       # mean of all values
      variance <- varEnv  # variance of all values
      nugget <- 0      # how much noise
      scale <- 10       # correlated on which scale
      alpha <- 2       # how strongly corrlated (covariance)
      f <- GaussRF(x=gridSize, y= gridSize, model=model, grid=TRUE,
                    param=c(mean, variance, nugget, scale, alpha))
      while((min(f) < minEnv)|(max(f) > maxEnv)) f <- ifelse(f < minEnv, abs(f), ifelse (f > maxEnv, 2*maxEnv-f, f))
   }
   if (ltype == "lr"){
      ###### random landscape
      f <- matrix(runif(gridside*gridside, minEnv, maxEnv), ncol=gridside)
   }
   if (ltype == "lg"){
      ######  gradient landscape
      my_land_tmp <- matrix(rep(seq(minEnv, maxEnv, length.out = gridside), each=gridside), ncol=gridside)
      f <- matrix(rnorm(gridside*gridside, 0, 2), ncol=gridside)
      f <- my_land_tmp + f
      while((min(f) < minEnv)|(max(f) > maxEnv)) f <- ifelse(f < minEnv, abs(f), ifelse (f > maxEnv, 2*maxEnv-f, f))
   }
   return(f)
}

#ex.
#  gridSize <- 0:(10 - 1) # x and y range of the grid
#  land_r <- createLand("lr")
#  land_r<-round(land_r)
#  image(gridSize, gridSize, land_r)

  
