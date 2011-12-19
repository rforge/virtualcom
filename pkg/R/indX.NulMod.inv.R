indX.NulMod.inv <- function(parameters, com_inv, dist.phy, dist.fun, invID, ...){
  # to test: parameters <- paramNull[1,]
 
 	null.mod <- parameters["null.mod"]
  	invdOUT <-  parameters["invaderOut"]
  
  	if (sum(com_inv[, invID], na.rm=TRUE) !=0 ){ 	# if there was at least one invasion
  	
  	
  	#======================================
  	# for random invasion null models
  	#======================================  	
  	
    	if (null.mod %in% c(1,2)){ 					
    	  	res <- lapply(invID, function(x){		# pour chaque invasive
     			tmp <- com_inv
				tmp[,x] <- 1							
     			tmp_inv <- div.param.invasion.obs(tmp, dist.phy, dist.fun, x)	# Compute the indices for all sites
     			tmp_sites <- as.data.frame(cbind(indiv=rep(x,length(com_inv[,x])),y=ifelse(com_inv[,x]==0,0,1),tmp_inv[[x]]))
     			return(tmp_sites)
   		  	})
	     	out_sites <- ldply(res,data.frame)	
    	 	#weight <- rep(1,length(com_inv[,x])) # to make the constrained randominvasion
			out_sites$fun_MDNC <- as.numeric(as.character(out_sites$fun_MDNC))
			out_sites$fun_MDNN <- as.numeric(as.character(out_sites$fun_MDNN))
			out_sites$fun_MDWNC <- as.numeric(as.character(out_sites$fun_MDWNC))
			out_sites$fun_MDMAS <- as.numeric(as.character(out_sites$fun_MDMAS))
			
      		out_sites$fun_MDNCsq <- as.numeric(as.character(out_sites$fun_MDNC))^2
	  		out_sites$fun_MDNNsq <- as.numeric(as.character(out_sites$fun_MDNN))^2
	   		out_sites$fun_MDWNCsq <- as.numeric(as.character(out_sites$fun_MDWNC))^2
	   		out_sites$fun_MDMASsq <- as.numeric(as.character(out_sites$fun_MDMAS))^2
         	out_sites$SiteID <- rep(rownames(com_inv),length(invID))
         		
     		mods <- sapply(c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS"),function(i){
     			glm_inv <- glmer(y~eval(parse(text=i))+eval(parse(text=paste(i,"sq",sep=""))) + (1|SiteID) + (1|indiv), data=out_sites ,family=binomial)
			   sum_glm <- fixef(glm_inv)
			   slop <-  as.numeric(sum_glm[2])
			   slopsq <-  as.numeric(sum_glm[3])
			interc <- as.numeric(sum_glm[1])
			 pval_slp <- summary(glm_inv)@coefs[2,4]
			   pval_slpsq <- summary(glm_inv)@coefs[3,4]
			   pval_interc <- summary(glm_inv)@coefs[1,4]
			   Rsq <- as.numeric(glm_inv@deviance["ML"])
			return(c(slop=slop, slopsq=slopsq, interc=interc, pval_slp=pval_slp, pval_slpsq=pval_slpsq, pval_interc=pval_interc, Rsq=Rsq))
      		})
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
    	 		# if(null.mod == 4){out_sites$weight <- ifelse(out_sites[,1]==1 & out_sites[,2]=="I" | out_sites[,1]==0 & out_sites[,2]=="N",1,0)}		# for constrained inv_nat
     			if(null.mod == 5){out_sites$weight <- ifelse(out_sites[,1]==1 & out_sites[,2]=="I" | out_sites[,1]==0 & out_sites[,2]=="I",1,0)}		# for non-successfulinv
 	     		sub_out <- out_sites[out_sites$weight==1,]
 	     		return(sub_out)
 		    })	
 	     
 	    	all_sites <- ldply(res, data.frame)
 	 	 	all_sites$fun_MDNCsq <- all_sites$fun_MDNC^2
	 	  	all_sites$fun_MDNNsq <- all_sites$fun_MDNN^2
	 	  	all_sites$fun_MDWNCsq <- all_sites$fun_MDWNC^2
	 	  	all_sites$fun_MDMASsq <- all_sites$fun_MDMAS^2
    
    		mods <- sapply(c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS"),function(i){
       			glm_inv <- glmer(PresAbs~eval(parse(text=i)) + eval(parse(text=paste(i,"sq",sep=""))) + (1|siteID),data=all_sites ,family=binomial)
				sum_glm <- fixef(glm_inv)
    			slop <-  as.numeric(sum_glm[2])
     			slopsq <-  as.numeric(sum_glm[3])
      			interc <- as.numeric(sum_glm[1])
				pval_slp <- summary(glm_inv)@coefs[2,4]
				pval_slpsq <- summary(glm_inv)@coefs[3,4]
				pval_interc <- summary(glm_inv)@coefs[1,4]
				Rsq <- as.numeric(glm_inv@deviance["ML"])
				return(c(slop=slop, slopsq=slopsq, interc=interc, pval_slp=pval_slp, pval_slpsq=pval_slpsq, pval_interc=pval_interc, Rsq=Rsq))
	      	})
    	}    
  	} else stop ("Please contact Mars, no UFO have been recorded")
	
	return(list(null.mod=null.mod, res=res, mods=mods)) 
}
         
