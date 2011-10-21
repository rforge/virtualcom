mntd_median <- function (samp, dis){
    N <- dim(samp)[1]			# How many communities do we have
    mntd <- numeric(N)			# create a vector of 0
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            diag(sample.dis) <- NA
            mntd[i] <- median(apply(sample.dis, 2, min, na.rm = TRUE))
        } else mntd[i] <- NA
    }
    mntd
}   
