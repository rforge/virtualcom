create.taxa<-function(tree,taxa.level){
#This function groups species according to their phylogenetic proximity. 
#The level of grouping depends on taxa.level. 
#If taxa.level is between 0 and 1, it will cut the tree at the given proportion of the tree (0 at the leafs, 1 at the root)
#If taxa.level is an integer higher than 1, it will generate taxa.level groups.  
#If any element is higher than 1, the function round will insure that taxa.level becomes an integer vector.
#If the cutting process gives every species its own category or put everyone in the same one, the function returns "NULL"

  if (any(taxa.level>1)){taxa.level<-round(taxa.level)}
  taxa.level<-taxa.level[!duplicated(taxa.level)]
  treeh<-as.hclust.phylo(tree)
  if (any(taxa.level<1)){
  d<-max(cophenetic(treeh))
  groups<-as.data.frame(cutree(treeh,h=taxa.level*d))
  }else{groups<-as.data.frame(cutree(treeh,k=taxa.level))}
  
  colnames(groups)<-as.character(taxa.level)
  taxa<-lapply(1:length(taxa.level),function(x){as.factor(groups[,x])})
  nlimit<-c(1,length(tree$tip.label))
  taxa2<-lapply(taxa,function(x){
    if(any(nlimit==nlevels(x))){y<-NULL
          }else{y<-x}
    return(y)})
  names(taxa2)<-as.character(taxa.level)
  return(taxa2)
}
