#native_inv: within site replace invader with another species not present in the site

NC_native_inv <- function(samp, dis, inva){
        # select invaded sites
        natives <- colnames(samp)[!colnames(samp) %in% inva]
        new.invader <- sample(natives, 1)
        new.samp <- samp
        new.samp[,inva] <- 0
        new.samp[, new.invader] <- samp[,inva]
        return(list(samp=new.samp, inv=new.invader))
}
