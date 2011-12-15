#aggregate communities on landscape within groups of pixels using a moving window
#com_land: table with env value, x:y coordinates and species (communities on landscape)
#cell: size of side of moving window

moving.window<-function(com_land, cell=1, step.length=1){
  #cell=2
  
  #landscape dimensions
  side<-length(com_land[com_land$X==1,"X"])
  new.com <- seq(1, (side - (cell- step.length)), by=step.length)
  side2 <- length(new.com)
  nsp<-ncol(com_land)-3
  
  #aggregated landscape structure
  agg_land <- as.data.frame(matrix(NA,nrow=(side2*side2),ncol=ncol(com_land),dimnames=list(c(1:(side2*side2)),names(com_land))))
  agg_land$X<-rep(1:side2,each=side2)
  agg_land$Y<-rep(1:side2,side2)
  
  for (i in 1:side2)  {
    for (j in 1:side2)  {
        agg_land[which(agg_land$X==i & agg_land$Y==j),c(4:(nsp+3))] <-  colSums(com_land[which(com_land$X %in% c(new.com[i]:(new.com[i]+(cell-1))) & com_land$Y %in% c(new.com[j]:(new.com[j]+(cell-1)))),4:(nsp+3)])
        agg_land[which(agg_land$X==i & agg_land$Y==j),"env"] <-  mean(com_land[which(com_land$X %in% c(new.com[i]:(new.com[i]+(cell-1))) & com_land$Y %in% c(new.com[j]:(new.com[j]+(cell-1)))),"env"])
        }
      }
  
  return(agg_land)
}
