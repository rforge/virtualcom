indX.NulMod.inv <- function(parameters, com_inv, com_nat, param.com_nat=50 , dist.phy, dist.fun, invID, ...){
  # to test: parameters = paramNull[1,]; com_inv = com_inv ; dist.phy = dist.phy ;  
  # dist.fun = dist.fun ; invID = invID ; com_nat = com_nat ; param.com_nat=50
 
 	null.mod <- parameters["null.mod"]
  	invdOUT <-  parameters["invaderOut"]
  
  	if (sum(com_inv[, invID], na.rm=TRUE) !=0 ){ 	# if there was at least one invasion
  	
  	
  	#======================================
  	# for random invasion null models
  	#======================================  	
  	
    	if (null.mod %in% c(1,2)){ 				# 1: random invasion in resistant sites, 2: random invasion in any site
    			
    	  	res <- lapply(invID, function(x){		# pour chaque invasive
    	  		if(null.mod == 1)  	tmp <- com_inv
    	  		if(null.mod == 2)	tmp <- rbind(com_inv, com_nat[sample( ceiling(param.com_nat*dim(com_nat)[1]/100 )),])	# param.com_nat = the proportion of non-invaded native comm added
				tmp[,x] <- 1							
     			tmp_inv <- div.param.invasion.obs(tmp, dist.phy, dist.fun, x)	# Compute the indices for all sites
     			tmp_sites <- as.data.frame(cbind(indiv=rep(x,length(com_inv[,x])),y=ifelse(com_inv[,x]==0,0,1),tmp_inv[[x]]))
     			return(tmp_sites)
   		  	})
	     		     	
	     	all_sites <- ldply(res,data.frame)	
    	 	#weight <- rep(1,length(com_inv[,x])) # to make the constrained randominvasion
			all_sites$y <- as.numeric(as.character(all_sites$y))
			
			all_sites$fun_MDNC <- as.numeric(as.character(all_sites$fun_MDNC))
			all_sites$fun_MDNN <- as.numeric(as.character(all_sites$fun_MDNN))
			all_sites$fun_MDWNC <- as.numeric(as.character(all_sites$fun_MDWNC))
			all_sites$fun_MDMAS <- as.numeric(as.character(all_sites$fun_MDMAS))
			all_sites$phy_MDNC <- as.numeric(as.character(all_sites$phy_MDNC))
			all_sites$phy_MDNN <- as.numeric(as.character(all_sites$phy_MDNN))
			all_sites$phy_MDWNC <- as.numeric(as.character(all_sites$phy_MDWNC))
			all_sites$phy_MDMAS <- as.numeric(as.character(all_sites$phy_MDMAS))
			
      		all_sites$fun_MDNCsq <- as.numeric(as.character(all_sites$fun_MDNC))^2
	  		all_sites$fun_MDNNsq <- as.numeric(as.character(all_sites$fun_MDNN))^2
	   		all_sites$fun_MDWNCsq <- as.numeric(as.character(all_sites$fun_MDWNC))^2
	   		all_sites$fun_MDMASsq <- as.numeric(as.character(all_sites$fun_MDMAS))^2
         	all_sites$SiteID <- rep(rownames(com_inv),length(invID))
        }
    
   	#======================================
  	# for swap invader null models
  	#======================================

    	if(null.mod %in% c(3:5)){	# 3: swap inva_nat, 4: constrained swap inv_nat 5: swap_inv_nonsuccessfulinv
    		tmp <- com_inv
    		res <-  lapply(c(1:dim(tmp)[1]), function(x){		# pour tous les sites
         		tabl <- sapply(colnames(tmp), function(y){	# MDNC pour toutes les especes
                	tmp <- com_inv
					tmp[,y] <- 1		
        	 		tmp_inv <- div.param.invasion.obs(rbind(tmp[x,],tmp[x,]), dist.phy, dist.fun, y)
         			return(tmp_inv)
    		 	})
    	 		names(tabl) <- sapply(strsplit(names(tabl),split=".",fixed=TRUE),function(x){x[1]})
    	 		tabl2 <- as.data.frame(t(sapply(tabl,function(x){x[1,]})))
    	 	
    	 		# Identify the sites that should be selected:
    	 		out_sites <- as.data.frame(cbind(ifelse(t(tmp[x,])==0,0,1),exoID=ifelse(row.names(tabl2)%in% invID,"I","N"),siteID=rep(x,dim(tabl2)[1]),tabl2))
     			names(out_sites)[1] <- "PresAbs"
     			if(null.mod == 3){out_sites$weight <- ifelse(out_sites[,1]==1 & out_sites[,2]=="I" | out_sites[,1]==0 & out_sites[,2]=="N",1,0)}		# for inv_nat
     			if(null.mod == 5){out_sites$weight <- ifelse(out_sites[,1]==1 & out_sites[,2]=="I" | out_sites[,1]==0 & out_sites[,2]=="I",1,0)}		# for non-successfulinv
     			sub_out <- out_sites[out_sites$weight==1,]
 	     		return(sub_out)
 		    })	
 	     
 	    	all_sites <- ldply(res, data.frame)
			all_sites$PresAbs <- as.numeric(as.character(all_sites$PresAbs))

 	    	all_sites$fun_MDNC <- as.numeric(as.character(all_sites$fun_MDNC))
			all_sites$fun_MDNN <- as.numeric(as.character(all_sites$fun_MDNN))
			all_sites$fun_MDWNC <- as.numeric(as.character(all_sites$fun_MDWNC))
			all_sites$fun_MDMAS <- as.numeric(as.character(all_sites$fun_MDMAS))

 	 	 	all_sites$fun_MDNCsq <- as.numeric(as.character(all_sites$fun_MDNC^2))
	 	  	all_sites$fun_MDNNsq <- as.numeric(as.character(all_sites$fun_MDNN^2))
	 	  	all_sites$fun_MDWNCsq <- as.numeric(as.character(all_sites$fun_MDWNC^2))
	 	  	all_sites$fun_MDMASsq <- as.numeric(as.character(all_sites$fun_MDMAS^2))
       	}    
  	} else stop ("Please contact Mars, no UFO have been recorded")
	
	return(list(null.mod = parameters, all_sites = all_sites)) 
}
         
