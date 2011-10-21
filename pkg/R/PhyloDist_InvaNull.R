################################################################################################
# R functions to compare the distance (phylogenetic, functional, ecological) between
# some species of interest (invasive, rare) and their co-occuring species in the communities with nulll models 
# and estimate how significant the deviation from random patterns is.
#
#
# INPUTS:
# - "samp" = Site x species matrix. A community matrix (abundance or presence/absence). Species names must be consistent with "dis". 
#
#
# - "dis" = A square distance matrix (make sure it is not asymetric). If it is just : as.matrix(dis). Comes from an ultrametric tree. 
#           Make sure you have the same species names (both for columns and rows) than in the "samp" matrix. 
#
#
# - "inva"= Vector of species names of interest. They are the species for which we'll estimate the distance with the communities
#			For instance, names of the invasive species in the pool of species. Make sure names match the names of the samp matrix.
#
# - "nrep"= Number of randomizations 
#
# - "abun" = Logical. Is there abundance data in the site x species matrix? If only presence/absence data only MDNC and MDNN are computed 
#
# - "nullModel" = which null model should be used to compare observed patterns. One of: 
#	"rand_inv": random invasions - shuffle sites for alien species - new previously non-invaded sites are invaded
#	"rand_invb": random invasions - shuffle sites for alien species - all sites may be invaded, including previously invaded sites
#	"native_inv": within site replace invader with a native not present in the site
#	"shuffle_dist": shuffle distance matrix/tree tips
#
# - "invaderOut" = Logical. should other invasive species be removed from the sample before calculating phylogenetic distance to the native community?
#		   Default is False (other invasive species are assumed to replace natives).
#
#
# OUTPUTS:
# list of 2 elements:
#
#	- datas: a list of lists in which each element is an invasive species. For each invasive species a list of dataframes for each metric. Dataframes have sites by rows
#		   and permutations by columns. Two rows with mean and median of the randomizations are added. 3 columns with the observed value, the frequency with which
#		   the observed > randomized, and the z-score are added.
#
#
#	- result: a list of summary tables for each metric. Invasive species by columns and summary stats by rows: number of invaded sites, p_value based on mean over sites,
#		    p_value based on median over sites, Number of sites with p-values below 0.05, number of sites with significant z-values, mean z and median z.
#
#
# METRICS CALCULATED:
#
#       - MDNN: Phylogenetic/functional/ecological distance of all species of interest to their nearest relative (native for instance).
#				This is done for each community in which the species of interest occur. 
#       - MDNC: Mean phylogenetic/Functional/ecological distance all species of interest to the native community.
#				This is done for each community in which the species of interest occur
#       - MDWNC: Weghted Mean phylogenetic/Functional/ecological distance all species of interest to the native community.
#				This is done for each community in which the species of interest occur. The mean distance is weighted by the abundance of each species
#       			The last two indices are somehow similar to the alpha component of Ackerly. 
#       - MDMAS: Mean distance to the most abundant species. If several species with the same maximum abundance, the mean is taken. 
#
#  
#
##########################################################################################


##########################################################################################
#  NULL MODELS



#native_inv: within site replace invader with another species not present in the site

native_inv <- function(samp, dis, inva, nrep=1, metric, invaderOut=FALSE){
    out <- list()
    temp.rand <- function(samp, dis, inva){
        # select invaded sites
        new.samp.tmp <- samp[eval(parse(text=paste("samp$", inva, ">0", sep=""))),]
        # change species by site matrix: 
	# For each invaded site set invader to 0 (absent) and another random species to 2 (later turned back to 1-present)
        new.samp <- t(apply(new.samp.tmp, 1, function(x) {
                                new.inv <-  sample(names(x[x == 0]), 1)
                                x[new.inv] <- 999
                                x[inva] <- 0
                                x
                                }))
        new.invader <- apply(new.samp, 1, function(x) names(x[x==999])) #keep track of the natives used as new invaders
        new.samp <- as.data.frame(ifelse(new.samp == 999, 1, new.samp))
        # to test: i=2; rbind(new.samp.tmp[i,],new.samp[i,])
        out.tmp <- list(vector())
        for (i in rownames(new.samp)){
           for(m in metric){
              tmp <- MDNC_MDNN_Invasive_n(new.samp[i,], dis, new.invader[i], invaderOut) # calculate in 1 site phylogenetic distance of 1 new invader
              out.tmp[[m]][i] <- tmp[[m]][1,1]
            }
        }
        return(out.tmp)
    } 
    for(x in inva){ # loop over species
       for(z in 1:nrep){ # loop over repetitions
          tmp.tmp <- temp.rand(samp, dis, x)
          for (y in metric){
             if(z==1) out[[x]][[y]] <- matrix(NA, ncol=nrep, nrow=nrow(samp), dimnames=list(rownames(samp),1:nrep))
             ind <- names(tmp.tmp[[y]])
             out[[x]][[y]][ind,z] <- tmp.tmp[[y]]
          }
          print(paste("species=",x,"; repetition=",z, sep=""))
       }
    }
    return(out)
}

#rand_inv: random invasions - shuffle sites for alien species - new previously non-invaded sites are invaded

rand_inv <- function(samp, dis, inva, nrep=1, metric, invaderOut=FALSE){
    # only non invaded sites selected (returns same values to sp that have not enough non invaded sites)
    out <- list()
	  # calculate for each invader MDNC in all sites and then use an index to randomly select invaded sites to take into account
    	  for (x in inva){
        	only.invader <- samp[,x]
        	tmp.samp <- cbind(1, samp[,!names(samp)==x]) 
       		names(tmp.samp)[1]=x
        	if (invaderOut==FALSE) tmp <- MDNC_MDNN_Invasive_n(tmp.samp, dis, x)
        	if (invaderOut==TRUE) tmp <- MDNC_MDNN_Invasive_n(tmp.samp, dis, inva,invaderOut) #if other aliens are to be removed need to pass list of species
	  	for (z in 1:nrep){
	            if(sum(only.invader) < sum(only.invader==0)) ind <- sample(length(only.invader), sum(only.invader), prob=ifelse(only.invader==1,0,1)) else  #sites already invaded get proability 0 of being resampled
	            ind <- which(only.invader==1) #if not enough sites, don't do randomization
	          	for (y in metric){
             		if(z==1) out[[x]][[y]] <- matrix(NA, ncol=nrep, nrow=nrow(samp), dimnames=list(rownames(samp),1:nrep))
             		out[[x]][[y]][ind, z] <- tmp[[y]][,x][ind]
         			}
     		    	print(paste("species=",x,"; repetition=",z, sep=""))
		    	}
    		}
    		return(out)
}


#rand_invb: random invasions - shuffle sites for alien species - all sites may be invaded, including previously invaded sites

rand_invb <- function(samp, dis, inva, nrep=1, metric, invaderOut=FALSE){
    # to test: samp=spxplot_2x2PA; dis= phyl.dist.nod.2x2; inva= neonat2x2; abun=FALSE;invaderOut=FALSE
    # x="Xanthium_orientale"
    out <- list()
    	  for (x in inva){
        		only.invader <- samp[,x]
        		tmp.samp <- cbind(1, samp[,!names(samp)==x])
       		names(tmp.samp)[1]=x
        		if (invaderOut==FALSE) tmp <- MDNC_MDNN_Invasive_n(tmp.samp, dis, x)
        		if (invaderOut==TRUE) tmp <- MDNC_MDNN_Invasive_n(tmp.samp, dis, inva)
  		for (z in 1:nrep){
            	ind <- sample(length(only.invader), sum(only.invader)) #all sites can be invaded
          		for (y in metric){
             		if(z==1) out[[x]][[y]] <- matrix(NA, ncol=nrep, nrow=nrow(samp), dimnames=list(rownames(samp),1:nrep))
             		out[[x]][[y]][ind, z] <- tmp[[y]][,x][ind]
         			}
     		    	print(paste("species=",x,"; repetition=",z, sep=""))
		    	}
    		}
    return(out)
}

# shuffle_dist: shuffle distance matrix/tree tips

shuffle_dist <- function(samp, dis, inva, nrep=1, metric, invaderOut=FALSE){
    out <- list()
  	for (z in 1:nrep){
		tmp.sp <- sample(names(dis))
    	  	tmp.dis <- dis
		colnames(tmp.dis)<-rownames(tmp.dis)<-tmp.sp
		tmp <- MDNC_MDNN_Invasive_n(samp, tmp.dis, inva, invaderOut=invaderOut)
		    for(x in inva){ # loop over species
       		for (y in metric){
             		if(z==1) out[[x]][[y]] <- matrix(NA, ncol=nrep, nrow=nrow(samp), dimnames=list(rownames(samp),1:nrep))
             		out[[x]][[y]][,z] <- tmp[[y]][,x]
         		 }
     		     }
		print(paste("repetition=",z, sep=""))
      }
    return(out)
}




################################################################################


# a function that calculates for each species the mean and median over sites of each randomization and of the observed pattern, for all invaded sites the frequency 
# with which the oberved pattern is greater than the randomized values (quantile_obs), Z value for each invaded site. Based on these values a result table is produced
# including number of invaded sites, p_value based on mean over sites, p_value based on median over sites, Number of sites with p-values below 0.05, number of sites with significant z-values, mean z and median z
# the function returns a list with the result table and for each species a list of dataframes for each metric with sites by rows and permutations by columns.

analyseRand <- function(inva, metric, datas){
# to test: datas=null
    result <- list()
    for (y in metric){
        result[[y]] <- matrix(NA, ncol=length(inva), nrow=9, dimnames=list(c("N_invaded","p_site_mean","p_site_median","N_sig_sites","N_sig_z","p_N_sig_sites","p_N_sig_z","mean_z","median_z"),inva))
        for (x in inva){
             colMean=colMeans(datas[[x]][[y]], na.rm=TRUE)
             colMedian=apply(datas[[x]][[y]], 2, median, na.rm=TRUE)
             #quantile_obs=apply(datas[[x]][[y]], 1, function(o) ifelse(sum(!is.na(o)) == 0, NA, sum(o[names(o)=="obs"] > o[!names(o)=="obs"], na.rm = TRUE) / length(o[!names(o)=="obs"])))
             quantile_obs=apply(datas[[x]][[y]], 1, function(o) ifelse(sum(!is.na(o)) == 0, NA, mean(o[names(o)=="obs"] > o[!names(o)=="obs"], na.rm = TRUE)))
             z_value=apply(datas[[x]][[y]], 1, function(o) o[names(o)=="obs"] - mean(o[!names(o)=="obs"], na.rm = TRUE) / sd(o[!names(o)=="obs"], na.rm = TRUE))
             
             result[[y]][["N_invaded", x]] <- sum(!is.na(datas[[x]][[y]][,"obs"]))  
             result[[y]][["p_site_mean", x]] <- 1 - mean(colMean["obs"] > colMean[names(colMean)!="obs"], na.rm = TRUE)
             result[[y]][["p_site_median", x]] <- 1- mean(colMedian["obs"] > colMedian[names(colMedian)!="obs"], na.rm = TRUE)
             result[[y]][["N_sig_sites", x]] <- sum(quantile_obs > 0.95, na.rm = TRUE)
             result[[y]][["N_sig_z", x]] <- sum(pnorm(z_value, 0, 1) > 0.95, na.rm = TRUE)  
             result[[y]][["p_N_sig_sites", x]] <- 1 - mean(quantile_obs > 0.95, na.rm = TRUE)
             result[[y]][["p_N_sig_z", x]] <- 1 - mean(pnorm(z_value, 0, 1) > 0.95, na.rm = TRUE)
             result[[y]][["mean_z", x]] <- mean(z_value, na.rm = TRUE)
             result[[y]][["median_z", x]] <- median(z_value, na.rm = TRUE)

	     # Add to sitexpermutations dataframes: site means and medians, frequency with which
	     # the observed > randomized (quantile_obs), and the z-score 
             datas[[x]][[y]] <- rbind(datas[[x]][[y]], colMean)
             datas[[x]][[y]] <- rbind(datas[[x]][[y]], colMedian)
             datas[[x]][[y]] <- cbind(datas[[x]][[y]], quantile_obs = c(quantile_obs, NA, NA))
             datas[[x]][[y]] <- cbind(datas[[x]][[y]], z_value = c(z_value, NA, NA))
                                           
       }
    }                                       
    return(list(datas=datas,result=result))
} 

################################################################################

#A wrapper for the results, which calculates the observed patterns, calls the appropriate null model for comparison, and returns the
#combined randomization tables for each species and summary of results

PhyloDist_InvaNull <- function(samp, dis, inva, nrep=1, abun=FALSE, nullModel=c("rand_inv","rand_invb","native_inv","shuffle_dist"), invaderOut=FALSE){
   # calculate observed values
   obs <- MDNC_MDNN_Invasive_n(samp=samp, dis=dis, inva=inva, invaderOut)
   if (abun == FALSE) metric <- c("MDNC","MDNN")
   if (abun == TRUE) metric <- c("MDNC","MDNN","MDWNC","MDMAS")

   # create randomizations
   null <- switch(nullModel,
                  "rand_inv" = rand_inv(samp, dis, inva, nrep=nrep, metric=metric,invaderOut=invaderOut),
		  "rand_invb" = rand_invb(samp, dis, inva, nrep=nrep, metric=metric,invaderOut=invaderOut), 
		  "native_inv" = native_inv(samp, dis, inva, nrep=nrep, metric=metric,invaderOut=invaderOut),
		  "shuffle_dist" = shuffle_dist(samp, dis, inva, nrep=nrep, metric=metric,invaderOut=invaderOut))

   # combine randomizations with observations
   for (x in inva){
      for (y in metric){
          null[[x]][[y]] <- cbind(null[[x]][[y]][rownames(obs[[y]]),], obs=obs[[y]][,x])
      }
   }

   # analyse results
   return(analyseRand(inva, metric=metric, null))   
}

             
################################################################################
