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
