################################################################################
# distribute communities on landscapes
# spxenvxcom: array of env values by species by communities
# land: landscape of environmental values (matrix)

land.com<-function(spxenvxcom, land){

  #landscape table
  side<-dim(land)[1]
  #community array
  ncom<-dim(spxenvxcom)[3]
  nsp<-dim(spxenvxcom)[2]-1
  nenv<-dim(spxenvxcom)[1]
  
  #landscape table
  colnames(land)<-paste("x",1:side, sep="")
  rownames(land)<-paste("y",1:side, sep="")
  coor_env<- cbind(stack(as.data.frame(land)), y=paste("y",1:side, sep=""))
  names(coor_env)= c("env","X","Y")
  coor_env[,2]<-as.numeric(substr(coor_env[,2],2,4))
  coor_env[,3]<-as.numeric(substr(coor_env[,3],2,4))

  #merge landscape table with com array
  com_land<-matrix(NA, nrow=nrow(coor_env),ncol=dim(spxenvxcom)[2],dimnames=list(1:nrow(coor_env),dimnames(spxenvxcom)[[2]]))
  for (i in 1:nrow(coor_env)){
     com_land[i,] <- spxenvxcom[spxenvxcom[,1,1]==coor_env[i,1],,sample(1:ncom,1)]
  }
  com_land<-data.frame(coor_env,com_land[,-1])
  
  return(com_land)
}  
  
  


