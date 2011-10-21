################################################################################################
# MDNC_MDNN_Invasive is a R function to estimate the distance (phylogenetic, functional, ecological) between
# some species of interest (invasive, rare) and their co-occuring species in the communities. 
#
#Author: Wilfried Thuiller, LECA
#Date: 4 February 2009 
#
# MDNC: Mean phylogenetic/Functional/ecological distance of one given (introduced, rare) species (sp) to the native community
# MDWNC: Mean weighted phylogenetic/Functional/ecological distance of one given (introduced, rare) species (sp) to the native community
# MDNN: Phylogenetic/functional/ecological distance of one given (introduced, rare) species (sp) to its nearest relative (native for instance)
# MDMAS: Phylogenetic/functional/ecological distance of one given (introduced, rare) species (sp) to the most abundance species 
# 
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
# - "invaderOut" = logical. If True removes invasive species (species that should not be accounted for) from the sample before calculating 
#			distance to the native community. Default is False (other invasive species are assumed to replace natives)
#
#
# OUTPUTS:
# list of 4 dataframes:
#       - MDNN: Phylogenetic/functional/ecological distance of all species of interest to their nearest relative (native for instance).
#				This is done for each community in which the species of interest occur. 
#       - MDNC: Mean phylogenetic/Functional/ecological distance all species of interest to the native community.
#				This is done for each community in which the species of interest occur
#       - MDWNC: Weghted Mean phylogenetic/Functional/ecological distance all species of interest to the native community.
#				This is done for each community in which the species of interest occur. The mean distance is weighted by the abundance of each species
#       The last two indices are somehow similar to the alpha component of Ackerly. 
#       - MDMAS: Mean distance to the most abundant species. If several species with the same maximum abundance, the mean is taken. 
#
# These two dataframes have sites by rows and species of interest names by columns. 
#
##########################################################################################



MDNC_MDNN_Invasive_n <- function (samp, dis, inva, invaderOut=FALSE) 
{

#########################################	
#Function which estimate Mean phylogenetic/Functional distance of one given introduced species (sp) to the native comm and 
#Phylogenetic/Functional distance of one given introduced species (sp) to its nearest native relative
	mdnc_fun <- function(samp_temp, dis, inva, sp, invaderOut){
		N <- dim(samp_temp)[1]	#number of plots
		mdnc_sp <- mdnn_sp <- mdwnc_sp <- mdmas_sp <- numeric(N)
		for (i in 1:N) {
			#remove species not in the give community
        		sppInSample <- names(samp_temp[i, samp_temp[i, ] > 0]) 
        		#remove invasive species (species that should not be accounted for)
        		if(invaderOut==TRUE) sppInSample <- c(sp, sppInSample[-which(sppInSample %in% intersect(sppInSample, inva))]) #CHANGED: introduced variable invaderOut           
        		if (length(sppInSample) > 1) {
				#extract the distance between species "sp" and the remaining community        	
  		        sample.dis <- dis[sppInSample[-which(sppInSample==sp)], sp]
          		#Distance to the nearest species
              mdnn_sp[i] <-min(sample.dis)
           		#Distance the unweighted mean community
              mdnc_sp[i] <- mean(sample.dis)                                                     
              #distance to the weighted mean community
              mdwnc_sp[i] <- w.mean(sample.dis, samp_temp[i,sppInSample[-which(sppInSample==sp)]])
              #distance to the most abundant species
              mdmas_sp[i] <- most.abdt(sample.dis, samp_temp[i,sppInSample[-which(sppInSample==sp)]])                              
        		}
        		else {
            		mdnn_sp[i] <-  mdnc_sp[i] <- mdwnc_sp[i] <- mdmas_sp[i] <- 0				# CHANGED: added [i] for mdwnc_sp and mdmas_sp as well
        		}
      	}
      res=list()
      res$mdnn_sp=mdnn_sp
      res$mdnc_sp=mdnc_sp
      res$mdwnc_sp=mdwnc_sp
      res$mdmas_sp=mdmas_sp
      res		
	}


#Estime la moyenne pondérée par l'abondance de chaque espèce dans la communauté	
w.mean<-function(data, abun){
  res<-sum(data*(abun/sum(abun)))
  return(res)
}

most.abdt<-function(data, abun){
  res<-mean(data[abun==max(abun)])
}
	
##########################################
	
	MDNC <- MDNN <- MDWNC <- MDMAS <- as.data.frame(matrix(NA, ncol=length(inva), nrow=nrow(samp), dimnames=list(rownames(samp), inva))) # CHANGED: default was 0 before

	for(j in 1:length(inva)){
		samp_temp <- samp[eval(parse(text=paste("samp$", inva[j], ">0", sep=""))),]
		if(nrow(samp_temp)==0) stop(paste("The species", inva[j], "does not appear in any community", sep=" "))
		res <- mdnc_fun(samp_temp, dis, inva, inva[j], invaderOut)
		MDNC[rownames(samp_temp),j] <- res$mdnc_sp
		MDNN[rownames(samp_temp),j] <- res$mdnn_sp
		MDWNC[rownames(samp_temp),j] <- res$mdwnc_sp
		MDMAS[rownames(samp_temp),j] <- res$mdmas_sp		
	}
	result=list()
	result$MDNC=MDNC
	result$MDNN=MDNN
	result$MDWNC=MDWNC
	result$MDMAS=MDMAS
	return(result)
}
