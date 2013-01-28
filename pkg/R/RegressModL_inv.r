regressMod <- function(finl){
		#	finl <- finL[[1]]

	null.mod <- finl[[1]]
	all_sites <- finl[[2]]
	
	#======================================
  	# for random invasion null models
  	#======================================  	
  	
    if (null.mod %in% c(1,2)){ 	
    		
		mods <- lapply(c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS","phy_MDNC", "phy_MDNN", "phy_MDWNC", "phy_MDMAS"),function(i){
			
			tmp_allSite <- all_sites
			tmp_allSite$aa <- eval(parse(text=paste("tmp_allSite$",i,sep="")))/100
     		if (i %in% c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS")){
				tmp_allSite$bb <- eval(parse(text=paste("tmp_allSite$",i,"sq",sep="")))/100
			} else { tmp_allSite$bb <- (eval(parse(text=paste("tmp_allSite$",i,sep="")))^2)/100 }
			
     		glm_inv12 <- glmer(y ~ aa + bb + (1|SiteID) + (1|indiv), data=tmp_allSite ,family="binomial")
     		glm_inv1 <- glmer(y ~ aa +  (1|SiteID) + (1|indiv), data=tmp_allSite ,family="binomial")
     		glm_inv2 <- glmer(y ~ bb +  (1|SiteID) + (1|indiv), data=tmp_allSite ,family="binomial")
     		glm_inv0 <- glmer(y~ 1 + (1|SiteID) + (1|indiv), data=tmp_allSite ,family="binomial")
     			
     		# par(mfrow=c(2,4))
     		# for( i in c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS","phy_MDNC", "phy_MDNN", "phy_MDWNC", "phy_MDMAS")){	vioplot(eval(parse(text=paste("all_sites[which(all_sites$y==0),'",i,"']",sep=""))),eval(parse(text=paste("all_sites[which(all_sites$y==1),'",i,"']",sep=""))),horizontal=TRUE,col="green")	}
     		# par(mfrow=c(1,2))
     		# plot(UnderInv~niche_evol,data=Pool$func,xlab="trait values",ylab="Potentially invase species")
     		# vioplot(Pool$func[Pool$func$UnderInv==0,"niche_evol"],Pool$func[Pool$func$UnderInv==1,"niche_evol"],breaks=100,col="white")	
     		
     		RES <- matrix(NA,3,7,dimnames=list(c("glm_inv1","glm_inv2","glm_inv12"),c("slop","slopsq","pval_slp","pval_slpsq","Radj","AUC","AIC")))
     		
     		for(j in c("glm_inv1","glm_inv2","glm_inv12"))	{
     			glm_inv  <- eval(parse(text=j))
     			sum_glm  <- fixef(glm_inv)
     			RES[j,"slop"] <- ifelse(j %in% c("glm_inv1","glm_inv12"),as.numeric(sum_glm[2]),NA)
     			RES[j,"slopsq"] <- ifelse(j=="glm_inv2",as.numeric(sum_glm[2]),ifelse(j=="glm_inv12",as.numeric(sum_glm[3]),NA))
     			RES[j,"pval_slp"] <- summary(glm_inv)@coefs[2,4]
     			RES[j,"pval_slpsq"] <- ifelse(dim(summary(glm_inv)@coefs)[1]>2,summary(glm_inv)@coefs[3,4],NA)
     			RES[j,"Radj"] <- 1-((summary(glm_inv)@logLik[[1]]-2)/summary(glm_inv0)@logLik[[1]])
     			RES[j,"AUC"] <- somers2(fitted(glm_inv),all_sites$y)["Dxy"]
     			RES[j,"AIC"] <- summary(glm_inv)@AICtab[1,1]
     		}	
			return(round(RES,5))
      	})
      	names(mods) <- c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS","phy_MDNC", "phy_MDNN", "phy_MDWNC", "phy_MDMAS")
    }
				

	#======================================
  	# for swap invader null models
  	#======================================

    if(null.mod %in% c(3:5)){	# 3: swap inva_nat, 4: constrained swap inv_nat 5: swap_inv_nonsuccessfulinv

		mods <- lapply(c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS","phy_MDNC", "phy_MDNN", "phy_MDWNC", "phy_MDMAS"),function(i){
						
			all_sites$aa <- eval(parse(text=paste("all_sites$",i,sep="")))/100
     		if (i %in% c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS")){
				all_sites$bb <- eval(parse(text=paste("all_sites$",i,"sq",sep="")))/100
			} else { all_sites$bb <- (eval(parse(text=paste("all_sites$",i,sep="")))^2)/100 }
			
			glm_inv12 <- glmer(PresAbs~ aa + bb + (1|siteID), data=all_sites, family="binomial")
			glm_inv1 <- glmer(PresAbs~ aa + (1|siteID), data=all_sites, family="binomial")
			glm_inv2 <- glmer(PresAbs~ bb + (1|siteID), data=all_sites, family="binomial")
       		glm_inv0 <- glmer(PresAbs~ 1 + (1|siteID), data=all_sites, family="binomial")

     		RES <- matrix(NA,3,7,dimnames=list(c("glm_inv1","glm_inv2","glm_inv12"),c("slop","slopsq","pval_slp","pval_slpsq","Radj","AUC","AIC")))
     		
     		for(j in c("glm_inv1","glm_inv2","glm_inv12"))	{
     			glm_inv  <- eval(parse(text=j))
     			sum_glm  <- fixef(glm_inv)
     			RES[j,"slop"] <- ifelse(j %in% c("glm_inv1","glm_inv12"),as.numeric(sum_glm[2]),NA)
     			RES[j,"slopsq"] <- ifelse(j=="glm_inv2",as.numeric(sum_glm[2]),ifelse(j=="glm_inv12",as.numeric(sum_glm[3]),NA))
     			RES[j,"pval_slp"] <- ifelse(j %in% c("glm_inv1","glm_inv12"), summary(glm_inv)@coefs[2,4],NA)
     			RES[j,"pval_slpsq"] <- ifelse(j=="glm_inv2",summary(glm_inv)@coefs[2,4],ifelse(j=="glm_inv12",summary(glm_inv)@coefs[3,4],NA))
     			RES[j,"Radj"] <- 1-((summary(glm_inv)@logLik[[1]]-2)/summary(glm_inv0)@logLik[[1]])
     			RES[j,"AUC"] <- somers2(fitted(glm_inv),all_sites$PresAbs)["Dxy"]
     			RES[j,"AIC"] <- summary(glm_inv)@AICtab[1,1]
     		}	
			return(round(RES,5))
		})
      	names(mods) <- c("fun_MDNC", "fun_MDNN", "fun_MDWNC", "fun_MDMAS","phy_MDNC", "phy_MDNN", "phy_MDWNC", "phy_MDMAS")
	}
	return(list(null.mod=null.mod, all_sites=all_sites, mods=mods)) 
}
    	     			











