constrained.NM<-function(spxp,taxa=NULL,sp.suit=NULL){
  spxp<-as.data.frame(spxp)
  if (ncol(spxp)==1){spxp<-t(spxp)}
  if (!is.null(taxa)){if (ncol(spxp)!=length(taxa)){stop("Non convenient taxa dimension ")}}
  if (!is.null(sp.suit)){
    if (any(dim(sp.suit)!=dim(spxp))){
      stop("Non convenient sp.suit dimensions")
    }
    sp.suit<-as.data.frame(sp.suit)
    if (ncol(sp.suit)==1){sp.suit<-t(sp.suit)}}

#Basic null model : tree tips shuffling
  if (is.null(taxa) & is.null(sp.suit)){
    spxp.rand<-spxp[,sample(1:ncol(spxp))]
    if(ncol(spxp.rand)==1){spxp.rand<-t(spxp.rand)}
  }

#Habitat scaled null model
  if (!is.null(sp.suit) & is.null(taxa)){
	site.rand<-function(site.suitability){
  		n<-length(site.suitability)/2
  		x<-site.suitability[1:n]
  		site.suit<-site.suitability[(n+1):(2*n)]
  		site.suit[site.suit==0]<-min(site.suit)*0.0001
  		site.compo<-rep(0,n)
  		values.pos<-x[x>0]
  		if (length(values.pos)>0){
  		sp.rand<-sample(1:n,length(values.pos),prob=site.suit)
  		site.compo[sp.rand]<-values.pos}
  	return(site.compo)
	}
      spxp.rand<-t(apply(cbind(spxp,sp.suit),1,site.rand))
  }

#Evolutionary scaled null model
  if (is.null(sp.suit) & !is.null(taxa)){
    n<-ncol(spxp)
    sp.split<-split(1:n,taxa)
    sp.rand<-unsplit(lapply(sp.split,function(x){
      if (length(x)>1){y<-sample(x,length(x))
                       }else{y<-x}
      return(y)
      }),taxa)
    spxp.rand<-as.data.frame(spxp[,sp.rand])
    if(ncol(spxp.rand)==1){spxp.rand<-t(spxp.rand)}
  }

#Habitat and Evolutionary scaled null model
  if (!is.null(sp.suit) & !is.null(taxa)){
	site.rand<-function(site.suitability){
  		n<-length(site.suitability)/2
  		x<-site.suitability[1:n]
  		site.suit<-site.suitability[(n+1):(2*n)]
  		site.suit[site.suit==0]<-min(site.suit)*0.0001
  		site.compo<-rep(0,n)
  		values.pos<-x[x>0]
  		if (length(values.pos)>0){
  		sp.rand<-sample(1:n,length(values.pos),prob=site.suit)
  		site.compo[sp.rand]<-values.pos}
  	return(site.compo)
	}
      tab.temp<-cbind(spxp,sp.suit)
      spxp.rand<-apply(tab.temp,1,function(y){
        y.split<-split(y,c(taxa,taxa))
        rand.split<-lapply(y.split,function(y.gr){site.rand(y.gr)})
        y.rand<-unsplit(rand.split,c(taxa,taxa))
        return(y.rand)})    
  }
  spxp
  rownames(spxp.rand)<-rownames(spxp)
  colnames(spxp.rand)<-colnames(spxp)  
  return(spxp.rand)
}
