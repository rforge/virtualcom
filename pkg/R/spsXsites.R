.spsXsites <- function(i, Nsp, my.matrix, my.list){
	for(j in 1:Nsp){ 		# To fill the siteXspecies matrix
 		if(colnames(my.matrix)[j] %in% colnames(my.list[[i]]$spcXs)) {
 			my.matrix[i,j] <- my.list[[i]]$spcXs[,colnames(my.matrix)[j]]
 		}                                                                                         
 	}
 	return(my.matrix[i,])
}
