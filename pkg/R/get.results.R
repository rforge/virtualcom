#' result collection
#' 
#' extracts results from the list that is given by \code{\link{simulation.experiment}}.
#' 
#' @param output list of results as it is given by simulation experiment
#' @param myvar variables that should be extracted. can be either "obs" for observed values, "zNULL" for z-values or "rankNULL" for ranks
#' @param invader should variables be extracted for the native or the invasive community? Default is FALSE and thus extraction for native community
#' @return data.frame
#' @seealso
#' \code{\link{simulation.experiment}} for running simulation experiments with an integrated simulation of species pools and diversity analyses of final community structures
#' 
#' @examples
#' # load pre-prepared parameter table
#' data(simple_param)	
#' # setting up a full experiment (this may take a few minutes) 
#' wrapper <- function(a){
#'   library(VirtualCom)
#'   data(simple_param)
#'   print(a)
#'   return(try(simulation.experiment(simple_param[a,])))	
#'   }  
#'   # running the experiment either with lapply
#'   output <- lapply(1:nrow(simple_param), wrapper) 
#'   ## running the experiment with sfLapply using snowfall to allow for parallel computing
#'   # require(snowfall)
#'   # sfInit(parallel=TRUE, cpus=11)
#'   # output <- sfLapply(1:nrow(simple_param), wrapper)
#'   # sfStop()
#'   
#'   # extract results from list
#'   result.table <- get.results(output=output, myvar="obs", invader="FALSE") 
#'   result.table$process <- ifelse(result.table$beta.env==0 & result.table$beta.comp==0, "Random", ifelse(result.table$beta.env!=0 & result.table$beta.comp==0, "Env", ifelse(result.table$beta.env==0 & result.table$beta.comp!=0, "Comp", "Both") ))
#'   
#'   # plot the funtional diversity (mpd) of communities in dependence on assembly processes
#'   require(ggplot2)
#'   ggplot(data=result.table, aes(x=process, y=FD_ab_mpd)) + geom_boxplot(width=0.8) +  xlab("Assembly rules") +  ylab("Functional diversity")
#'   
#' @export
 
get.results <- function(output, myvar, invader = "FALSE") {
    if (invader == "FALSE") {
        list.of.results <- lapply(1:length(output), function(x) {
            data.frame(c(output[[x]]$parameter, output[[x]]$pool$indices), eval(parse(text = paste("output[[", x, "]]$natives$indices$", myvar, sep = ""))))
        })
    } else {
        list.of.results <- lapply(1:length(output), function(x) {
            allInv <- names(output[[x]]$invaders$indices)
            if (length(allInv) >= 1) {
                for (inv in 1:length(allInv)) {
                  temp <- data.frame(c(output[[x]]$parameter, output[[x]]$pool$indices), eval(parse(text = paste("output[[", x, "]]$invaders$indices[[", inv, "]]$", myvar, 
                    sep = ""))), allInv[inv])
                  if (inv == 1) 
                    tempA <- temp else tempA <- rbind(tempA, temp)
                }
                if (nrow(tempA) >= 2) {
                  return(tempA)
                }
            }
        })
    }
    return(do.call("rbind", list.of.results))
} 
