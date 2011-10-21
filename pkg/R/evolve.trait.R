###############################################################################
# function for trait evolution (from R-Forge; Ackerly)
###############################################################################
# This function explicitly simulates the steps from one node to the other; i.e. it devides the branch between
# ancestor and daugther nodes into substeps and explicitly simulates changes from substep to substeps
# Parameter choices for different models ('trait.mode'):

# always:
#       x.root = 0          #root value
#       set.minbl = 1,      #length of shortest branch
#     	show.plots = TRUE   #plotting?
#     	use.color = FALSE   #plots with color
#     	debug = FALSE       #debugging?


# for 'Brownian':
#     	sigma = 1           #brownian motion st. dev.
#     	trend = 0           #brownian motion trend (gives the mean of the random variable)
#       burst = rep(1,nrow(phy$edge)) # edge-specific rate constant (constant development in one direction)


# for 'Bounded' (bounded Brownian motion):
#     	sigma = 1            #brownian motion st. dev.
#     	trend = 0            #brownian motion trend (gives the mean of the random variable)
#       bound = c(-10,10)    #bounds on trait evolution


# for 'OU' (ornstein-uhlenbeck)
#     	sigma = 1            #brownian motion st. dev.
#     	trend = 0            #brownian motion trend (gives the mean of the random variable)
#     	mu = rep(0,nrow(phy$edge))   #ornstein-uhlenbeck mode (attractor for each node)
#     	theta = 1                    #ornstein-uhlenbeck parameter (strength of attraction)


# for 'MH-OU' (mh-ornstein-uhlenbeck; cf. mu)
#     	sigma = 1            #brownian motion st. dev.
#     	trend = 0            #brownian motion trend (gives the mean of the random variable)
#     	mu = rep(0,nrow(phy$edge))   #ornstein-uhlenbeck mode (attractor for daugther nodes is the value of ancestor)
#     	theta = 1                    #ornstein-uhlenbeck parameter (strength of attraction)

# for 'ACDC'
#     	sigma = 1            #sigma changes for each substep and over time (new sigma = sigma/gamma)
#     	trend = 0            #brownian motion trend (gives the mean of the random variable)
#     	gamma = 1            #age shift factor for ACDC
#     	gamma.reset = 0      #probability of resetting sigma under ACDC


# for 'Proportional' (changes not additive but multiplicative with lognormal instead of normal noise)
#     	sigma = 1           #log normal distribution st. dev.
#     	trend = 0           #log normal distribution trend (gives the mean of the random variable)


# for 'Speciational' (only one speciation event in one substep; branch length not important)
#     	sigma = 1           #brownian motion st. dev.
#     	trend = 0           #brownian motion trend (gives the mean of the random variable)


# for 'Punctuational' (only one speciation event in one substep; only one daugther node with speciation event; branch length not important)
#     	sigma = 1           #brownian motion st. dev.
#     	trend = 0           #brownian motion trend (gives the mean of the random variable)

`evolve.trait` <-
function(phy,
	x.root=0, #root value
	sigma = 1, #brownian motion st. dev.
	trend = 0, #brownian motion trend
	bound = c(-10,10), # bounds on trait evolution
	burst = rep(1,nrow(phy$edge)), # edge-specific rate constant
	mu = rep(0,nrow(phy$edge)), #ornstein-uhlenbeck mode
	theta = 1, #ornstein-uhlenbeck parameter
	gamma = 1, #age shift factor for ACDC
	gamma.reset = 0, #probability of resetting sigma under ACDC
	trait.mode = 'Brownian',
	set.minbl = 1, #length of shortest branch
	show.plots.tree = FALSE,
	show.plots.trait = FALSE,
	use.color = FALSE, # plots with color
	debug = FALSE, ...) {
  # phy=myTree; trait.mode="OU"; x.root=0; sigma=sigma; mu=mymu; theta=0; show.plots.trait = TRUE; use.color=TRUE; main="alpha=0";set.minbl = 1
  
	# I think this code requires the tree to be in cladewise order -- pdc
	phy <- reorder(phy, 'cladewise') # reorder.phylo {ape}
	#print(mu)
	phy$edge.length.original = phy$edge.length
	minbl = min(phy$edge.length)
	phy$edge.length = round(set.minbl*phy$edge.length/minbl)

	phy$Nterm = length(phy$tip.label) # number leaves
	phy$Nedge = nrow(phy$edge) # number of all branches
	phy$node.traits = rep(NA,(phy$Nterm+phy$Nnode))
	if (length(grep(trait.mode,'Proportional'))==1) x.root = 1
	phy$node.traits[phy$edge[1,1]] = x.root
	phy$ages = node.age(phy)$ages    # Phylo object with phylo\$ages vector of node ages corresponding to phylo\$edge (cumulative!)
	max.ht = max(vcv.phylo(phy)) # This function computes the expected variances and covariances of a continuous trait assuming it evolves under a given model.
	min.trait = 0
	max.trait = 0

	if (length(grep(trait.mode,'ACDC'))==1) sg = rep(sigma,phy$Nedge)

	# Variable pu = 0 for first daughter of each node and = 1
	# for second daughter; used to assign character change to
	# one daughter in Punctuational mode
	if (length(grep(trait.mode,'Punctuational'))==1) {
		pu = rep(c(0,1),phy$Nnode)
		pu = pu[rank(phy$edge[,1],ties.method='first')]
	}

################################################################################
# here it starts:
	edges = list(NULL)
	for (i in 1:phy$Nedge) {
		edges[[i]] = matrix(NA, phy$edge.length[i]+1, 2)
		if (phy$edge[i,1]==phy$Nterm+1) {   # if it is root
			anc.age = 0
			anc.trait = x.root
			} else {
			anc.edge = match(phy$edge[i,1], phy$edge[,2]) # returns a vector of the positions of (first) matches of its first argument in its second
			anc.age = phy$ages[match(phy$edge[i,1],phy$edge[,2])] # gives cumulative age
			anc.trait = phy$node.traits[phy$edge[i,1]] # gives trait of ancestor
		}
		edges[[i]][,1] = anc.age:phy$ages[i] # devides time distance into steps?
		edges[[i]][1,2] = anc.trait
		if (length(grep(trait.mode,'Brownian'))==1) { # if Brownian motion
				  edges[[i]][-1,2] =
					anc.trait +
					cumsum(burst[i]*rnorm(phy$edge.length[i],
					mean = trend,
					sd = sigma))
			} else if (length(grep(trait.mode,'Bounded'))==1) { # if bounded Brownian motion
				for (j in 2:nrow(edges[[i]])) {
					repeat {
						edges[[i]][j,2] = edges[[i]][(j-1),2] +
							rnorm(1,mean=trend,sd=sigma)
						if ((edges[[i]][j,2] > bound[1]) &
							(edges[[i]][j,2] < bound[2])) break
					}
				}
			} else if (length(grep(trait.mode,'OU'))==1) {  # if OU process
				for (j in 2:nrow(edges[[i]])) {
					  #edges[[i]][j,2] = mu[i] +   #TODO: why?
						#theta*(edges[[i]][(j-1),2]-mu[i]) +
						#rnorm(1,mean=trend,sd=sigma)
						# I think this would be right:
					   edges[[i]][j,2] = edges[[i]][(j-1),2] +
						 theta*(mu[i] - edges[[i]][(j-1),2]) + 
             rnorm(1,mean=trend,sd=sigma)
				}
			} else if (length(grep(trait.mode,'MH-OU'))==1) { # if MH-OU process
				for (j in 2:nrow(edges[[i]])) {
					  edges[[i]][j,2] = mu[i] +   #TODO: why?
						theta*(edges[[i]][(j-1),2]-mu[i]) +
						rnorm(1,mean=trend,sd=sigma)
				}
				nextbranch = which(phy$edge[,1]==phy$edge[i,2])
				mu[nextbranch] = edges[[i]][nrow(edges[[i]]),2]
			} else if (length(grep(trait.mode,'ACDC'))==1) {# if ACDC
				sd.edge = sg[i] # sigma
				for (j in 2:nrow(edges[[i]])) {
					#varadj = gamma^(-edges[[i]][j,1])
					edges[[i]][j,2] = edges[[i]][(j-1),2] +
						rnorm(1,
						mean=trend,
						sd=sd.edge)
					sd.edge = sd.edge/gamma
					if (runif(1) < gamma.reset) sd.edge = sigma
				}
				nextbranch = which(phy$edge[,1]==phy$edge[i,2])
				sg[nextbranch] = sd.edge
			} else if (length(grep(trait.mode,	          # if Proportional
				'Proportional'))==1) {
				edges[[i]][-1,2] =
					anc.trait *
					cumprod(rlnorm(
					phy$edge.length[i],
					meanlog = trend,
					sdlog = sigma))
			} else if (length(grep(trait.mode,	          # if Speciational
				'Speciational'))==1) {
				edges[[i]][2,2] =
					anc.trait + rnorm(1,mean = trend,sd = sigma)
				edges[[i]][c(-1,-2),2] = edges[[i]][2,2]
			} else if (length(grep(trait.mode,
				'Punctuational'))==1) {                     # if Punctuational
				edges[[i]][2,2] =
					anc.trait + pu[i] *
					rnorm(1,mean = trend,sd = sigma)
				edges[[i]][c(-1,-2),2] = edges[[i]][2,2]
			} else {
				print('invalid trait evolution mode')
				break
			}
		if (min(edges[[i]][,2]) < min.trait) min.trait = min(edges[[i]][,2])
		if (max(edges[[i]][,2]) > max.trait) max.trait = max(edges[[i]][,2])

		phy$node.traits[phy$edge[i,2]] = edges[[i]][nrow(edges[[i]]),2]

		#edges = list(edges,edge)
	}
################################################################################

	x=matrix(phy$node.traits[1:phy$Nterm],ncol=1)
	row.names(x) = phy$tip.label

	if (show.plots.tree || show.plots.trait) {
#  		par(mfcol=c(2,1),cex.lab=1,cex.axis=1,
#  			mar=c(5.1,5.1,1.1,1.1))
  		phy.coph = cophenetic.phylo(phy)
  		trait.dist = as.matrix(dist(x,diag=TRUE,upper=TRUE))
  		# red=0, yellow=1/6, green=2/6, cyan=3/6, blue=4/6 and magenta=5/6
  		edgecolors <- rainbow(phy$Nedge,  start = 3/6, end = 5/6)
  		#edgecolors <- topo.colors(phy$Nedge)
  }
  if (show.plots.tree){
      #plot.phy <- phy
      #plot.phy$node.traits <- phy$node.traits
      if(use.color) {
          plot.phylo(phy, edge.color = edgecolors, show.tip.label=FALSE)
      } else {
          plot.phylo(phy, show.tip.label=FALSE)
      }
      ## TODO scale bar often overlaps with bottom branch
      # add.scale.bar(0,1,set.minbl)
  }
  if (show.plots.trait){
      plot(c(0,max.ht),c(min.trait,max.trait),type='n',xlab='',ylab='',...)
		  if(use.color) {
 			    for (i in 1:phy$Nedge) lines(edges[[i]], col = edgecolors[i])
		  } else {
 			    for (i in 1:phy$Nedge) lines(edges[[i]])
  		}
  }
	#print(phy$node.traits)
	#print(phy$edge)
	x <- as.vector(x[,1])
	names(x) <- phy$tip.label
	return(list(x, edges, c(min.trait,max.trait)))
}
