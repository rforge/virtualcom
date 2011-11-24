#aggregate communities on landscape within groups of pixels using a moving window
#com_land: table with env value, x:y coordinates and species (communities on landscape)
#cell: size of side of moving window

moving.window<-function(com_land, cell=1, step.length=1){
  #cell=2
  
  #landscape dimensions
  side<-length(com_land[com_land$X==1,"X"])
  new.com <- seq(1, (side - step.length), by=step.length)
  side2 <- length(new.com)
  
  #aggregated landscape structure
  agg_land <- as.data.frame(matrix(NA,nrow=(side2*side2),ncol=ncol(com_land),dimnames=list(c(1:(side2*side2)),names(com_land))))
  agg_land$X<-rep(1:side2,each=side2)
  agg_land$Y<-rep(1:side2,side2)
  
  for (i in 1:side2)  {
    for (j in 1:side2)  {
        agg_land[which(agg_land$X==i & agg_land$Y==j),c(4:103)] <-  colSums(com_land[which(com_land$X %in% c(new.com[i]:(new.com[i]+(cell-1))) & com_land$Y %in% c(new.com[j]:(new.com[j]+(cell-1)))),4:103])
        agg_land[which(agg_land$X==i & agg_land$Y==j),"env"] <-  mean(com_land[which(com_land$X %in% c(new.com[i]:(new.com[i]+(cell-1))) & com_land$Y %in% c(new.com[j]:(new.com[j]+(cell-1)))),"env"])
        }
      }
  
  return(agg_land)
}

##test
##single lands
#agg_land2<-moving.window(com_land_r, cell=2)
#agg_land3<-moving.window(com_land_r, cell=3)
#dim(agg_land3)
#agg_land2[1:2,1:10]
#agg_land3[1:2,1:10]
#
##all aggregated lands
#side<-length(com_land[com_land$X==1,"X"])
#agg_lands<-list()
#for(i in 1:side) agg_lands[[i]]<- moving.window(com_land_r, cell=i) 
#
##full landscape
#rbind(agg_lands[[10]][-2:-3],c(mean(com_land_r$env),colSums(com_land_r[4:103])))