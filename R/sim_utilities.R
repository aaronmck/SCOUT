#############################################################

#Citation: Zhang, X., Xu, C. & Yosef, N. Simulating multiple faceted variability in single cell RNA sequencing. Nat Commun 10, 2611 (2019). https://doi.org/10.1038/s41467-019-10500-w
#https://github.com/YosefLab/SymSim/blob/master/R/simulation_functions.R

#############################################################
# Master Equation Related Functions
#############################################################

#' Getting GeneEffects matrices 
#'
#' This function randomly generates the effect size of each evf on the dynamic expression parameters
#' @param ngenes number of genes
#' @param nevf number of evfs
#' @param randomseed (should produce same result if ngenes, nevf and randseed are all the same)
#' @param prob the probability that the effect size is not 0
#' @param geffect_mean the mean of the normal distribution where the non-zero effect sizes are dropped from 
#' @param geffect_sd the standard deviation of the normal distribution where the non-zero effect sizes are dropped from 
#' @param evf_res the EVFs generated for cells
#' @param is_in_module a vector of length ngenes. 0 means the gene is not in any gene module. In case of non-zero values, genes with the same value in this vector are in the same module.
#' @return a list of 3 matrices, each of dimension ngenes * nevf
GeneEffects <- function(ngenes,nevf,randseed,prob,geffect_mean,geffect_sd, evf_res, is_in_module){
  set.seed(randseed)
  
  gene_effects <- lapply(c('kon','koff','s'),function(param){
    effect <- lapply(c(1:ngenes),function(i){
      nonzero <- sample(size=nevf,x=c(0,1),prob=c((1-prob),prob),replace=T)
      nonzero[nonzero!=0]=rnorm(sum(nonzero),mean=geffect_mean,sd=geffect_sd)
      return(nonzero)
    })
    return(do.call(rbind,effect))
  })
  
  
  mod_strength <- 0.8
  if (sum(is_in_module) > 0){ # some genes need to be assigned as module genes; their gene_effects for s will be updated accordingly
    for (pop in unique(is_in_module[is_in_module>0])){ # for every population, find the 
      #which(evf_res[[2]]$pop == pop)
      for (iparam in c(1,3)){
        nonzero_pos <- order(colMeans(evf_res[[1]][[iparam]][which(evf_res[[2]]$pop==pop),])- colMeans(evf_res[[1]][[iparam]]),
                             decreasing = T)[1:ceiling(nevf*prob)]
        geffects_centric <- rnorm(length(nonzero_pos), mean = 0.5, sd=0.5)
        for (igene in which(is_in_module==pop)){
          for (ipos in 1:length(nonzero_pos)){
            gene_effects[[iparam]][igene,nonzero_pos[ipos]] <- rnorm(1, mean=geffects_centric[ipos], sd=(1-mod_strength))
          }
          gene_effects[[iparam]][igene,setdiff((1:nevf),nonzero_pos)] <- 0
        }
      }
      
      nonzero_pos <- order(colMeans(evf_res[[1]][[2]][which(evf_res[[2]]$pop==pop),])- colMeans(evf_res[[1]][[2]]),
                           decreasing = T)[1:ceiling(nevf*prob)]
      geffects_centric <- rnorm(length(nonzero_pos), mean = -0.5, sd=0.5)
      for (igene in which(is_in_module==pop)){
        for (ipos in 1:length(nonzero_pos)){
          gene_effects[[2]][igene,nonzero_pos[ipos]] <- rnorm(1, mean=geffects_centric[ipos], sd=(1-mod_strength))
        }
        gene_effects[[2]][igene,setdiff((1:nevf),nonzero_pos)] <- 0
      }
    }
  }
  return(gene_effects)
}

#' sample from smoothed density function
#' @param nsample number of samples needed
#' @param den_fun density function estimated from density() from R default
SampleDen <- function(nsample,den_fun){
  probs <- den_fun$y/sum(den_fun$y)
  bw <- den_fun$x[2]-den_fun$x[1]
  bin_id <- sample(size=nsample,x=c(1:length(probs)),prob=probs,replace=T)
  counts <- table(bin_id)
  sampled_bins <- as.numeric(names(counts))
  samples <- lapply(c(1:length(counts)),function(j){
    runif(n=counts[j],min=(den_fun$x[sampled_bins[j]]-0.5*bw),max=(den_fun$x[sampled_bins[j]]+0.5*bw))
  })
  samples <- do.call(c,samples)
  return(samples)
}

#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param match_param_den the fitted parameter distribution density to sample from 
#' @param bimod the bimodality constant
#' @param scale_s a factor to scale the s parameter, which is used to tune the size of the actual cell (small cells have less number of transcripts in total)
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params <- function(gene_effects,evf,match_param_den,bimod,scale_s){
  params <- lapply(1:3, function(iparam){evf[[iparam]] %*% t(gene_effects[[iparam]])})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    # X=matrix(data=c(1:10),ncol=2) 
    # this line is to check that the row and columns did not flip
    temp <- plyr::alply(X, 1, function(Y){Y})
    values <- do.call(c,temp)
    ranks <- rank(values)
    sorted <- sort(SampleDen(nsample=max(ranks),den_fun=match_param_den[[i]]))
    temp3 <- matrix(data=sorted[ranks],ncol=length(X[1,]),byrow=T)
    return(temp3)
  })
  
  bimod_perc <- 1
  ngenes <- dim(scaled_params[[1]])[2]; bimod_vec <- numeric(ngenes)
  bimod_vec[1:ceiling(ngenes*bimod_perc)] <- bimod
  bimod_vec <- c(rep(bimod, ngenes/2), rep(0, ngenes/2))
  scaled_params[[1]] <- apply(t(scaled_params[[1]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[2]] <- apply(t(scaled_params[[2]]),2,function(x){x <- 10^(x - bimod_vec)})
  scaled_params[[3]] <- t(apply(scaled_params[[3]],2,function(x){x<-10^x}))*scale_s
  
  return(scaled_params)
}


#' Getting the parameters for simulating gene expression from EVf and gene effects
#'
#' This function takes gene_effect and EVF, take their dot product and scale the product to the correct range 
#' by using first a logistic function and then adding/dividing by constant to their correct range
#' @param evf a vector of length nevf, the cell specific extrinsic variation factor
#' @param gene_effects a list of three matrices (generated using the GeneEffects function), 
#' each corresponding to one kinetic parameter. Each matrix has nevf columns, and ngenes rows. 
#' @param param_realdata the fitted parameter distribution to sample from 
#' @param bimod the bimodality constant
#' @return params a matrix of ngenes * 3
#' @examples 
#' Get_params()
Get_params2 <- function(gene_effects,evf,bimod,ranges){
  params <- lapply(gene_effects,function(X){evf %*% t(X)})
  scaled_params <- lapply(c(1:3),function(i){
    X <- params[[i]]
    temp <- apply(X,2,function(x){1/(1+exp(-x))})
    temp2 <- temp*(ranges[[i]][2]-ranges[[i]][1])+ranges[[i]][1]
    return(temp2)
  })
  scaled_params[[1]]<-apply(scaled_params[[1]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[2]]<-apply(scaled_params[[2]],2,function(x){x <- 10^(x - bimod)})
  scaled_params[[3]]<-apply(scaled_params[[3]],2,function(x){x<-abs(x)})
  scaled_params <- lapply(scaled_params,t)
  return(scaled_params)
}

get_prob <- function(glength){
  if (glength >= 1000){prob <- 0.7} else{
    if (glength >= 100 & glength < 1000){prob <- 0.78}
    else if (glength < 100) {prob <- 0}
  }
  return(prob)
}