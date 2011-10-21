div.param.invasion.obs <- function(ispSite, iphy, ifun, imyinva){
     phyloRes <- MDNC_MDNN_Invasive_n(samp=ispSite, dis=iphy, inva=imyinva)
     funRes <- MDNC_MDNN_Invasive_n(samp=ispSite, dis=ifun, inva=imyinva)

     iobs <- lapply(imyinva, function(x) {
     			tm <- cbind(sapply(names(phyloRes), function(y)  unlist(phyloRes[[y]][x])), sapply(names(funRes), function(y)  unlist(funRes[[y]][x]))); 
              	rownames(tm) <- 1: nrow(ispSite); 
              	colnames(tm) <- c(paste("phy_",names(phyloRes),sep=""),paste("fun_",names(phyloRes),sep=""));
              	tm})
     names(iobs) <- imyinva
     iobs
}

