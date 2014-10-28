get.attractors <- function(phy, theta, ngroups=2, method="single", plotting=FALSE, randomize=TRUE)
{
    if (ngroups > 1 & (!is.vector(theta) | length(theta) != ngroups)) print("theta has wrong fromat")
    # cluster tips and assign theta to tips; store everything in attra
    myDist <- as.dist(as.data.frame(cophenetic.phylo(phy)),upper=TRUE,diag=TRUE) # calculate distance matrix for tree
    myCluster <- hclust(myDist, method=method) #plot(myCluster)
    myGroups <- cutree(myCluster, k=ngroups)
    attra <- as.data.frame(cbind(node.index=match(names(myGroups), phy$tip.label), group=myGroups))
    if(randomize==TRUE) theta <- sample(theta)
    attra$attractor <- theta[attra$group]
    if(plotting==TRUE){
       phyTMP <- phy
       phyTMP$tip.label[attra$node.index] <- round(attra$attractor,1)
       parOld <- par(mfrow=c(1,2))
          plot(phyTMP)
          plot(phy)
       par(parOld)
    }

    # calculate theta for all inner nodes
    inner.nodes <- phy$edge[,2][phy$edge[,2] %in% phy$edge[,1]] # all inner nodes but root
    for (myNode in sort(inner.nodes, decreasing=TRUE)){  # over all inner nodes but root, start with outer nodes
         descendants <- phy$edge[,2][phy$edge[,1] == myNode] # descendants from node in focus
         myValue <- sample(attra$attractor[attra$node.index %in% descendants], size=1) # node in focus gets random value from descendants
         attra <- rbind(attra, c(myNode, NA, myValue))
    }
    # calculate theta for root
    descendantsRoot <- phy$edge[,2][phy$edge[,1] == min(phy$edge[,1])]
    attra <- rbind(attra, c(min(phy$edge[,1]), NA, sample(attra$attractor[attra$node.index %in% descendantsRoot]), size=1))
    return(attra)
}
