  library(tidyverse)
library(caret)
library(glmnet)
library(lmerTest)
library(igraph)
library(RcppAlgos)
library( prodlim)
library(doParallel)
library(tibble)
library("Rlab")  
library( randomForest)
library( dplyr)
library( partykit)
library( matrixStats)
library( grr)

.start_subset <- function(data) {
  ret <- 1:NROW(model.frame(data))
  if (length(data$yxmissings) > 0)
    ret <- ret[!(ret %in% data$yxmissings)]
  ret
}
constparties <- function(nodes, data, weights, fitted = NULL, terms = NULL, info = NULL) {
  
  stopifnot(all(sapply(nodes, function(x) inherits(x, "partynode"))))
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(weights, "list"))
  
  if(!is.null(fitted)) {
    stopifnot(inherits(fitted, "data.frame"))
    stopifnot(nrow(data) == 0L | nrow(data) == nrow(fitted))
    if (nrow(data) == 0L)
      stopifnot("(response)" %in% names(fitted))
  } else {
    stopifnot(nrow(data) > 0L)
    stopifnot(!is.null(terms))
    fitted <- data.frame("(response)" = model.response(model.frame(terms, data = data, 
                                                                   na.action = na.pass)),
                         check.names = FALSE)
  }
  
  ret <- list(nodes = nodes, data = data, weights = weights, fitted = fitted)
  class(ret) <- c("constparties", "parties")
  
  if(!is.null(terms)) {
    stopifnot(inherits(terms, "terms"))
    ret$terms <- terms
  }
  
  if (!is.null(info))
    ret$info <- info
  
  ret
}

.y2infl <- function(data, response, ytrafo = NULL) {
  
  if (length(response) == 1) {
    if (!is.null(ytrafo[[response]])) {
      yfun <- ytrafo[[response]]
      rtype <- "user-defined"
    } else {
      rtype <- class(data[[response]])[1]
      if (rtype == "integer") rtype <- "numeric"
      if (rtype == "AsIs") rtype <- "numeric"
    }
    response <- data[[response]]
    
    infl <- switch(rtype,
                   "user-defined" = yfun(response),
                   "factor" = {
                     X <- model.matrix(~ response - 1)
                     if (nlevels(response) > 2) return(X)
                     return(X[,-1, drop = FALSE])
                   },
                   "ordered" = {
                     sc <- attr(response, "scores")
                     if (is.null(sc)) sc <- 1L:nlevels(response)
                     sc <- as.numeric(sc)
                     return(matrix(sc[as.integer(response)], ncol = 1))
                   },
                   "numeric" = response,
                   "Surv" = .logrank_trafo(response)
    )
  } else {
    ### multivariate response
    infl <- lapply(response, .y2infl, data = data)
    tmp <- do.call("cbind", infl)
    attr(tmp, "assign") <- rep(1L:length(infl), sapply(infl, NCOL))
    infl <- tmp
  }
  if (!is.matrix(infl)) infl <- matrix(infl, ncol = 1)
  storage.mode(infl) <- "double"
  return(infl)
}

cforest_new <- function
(
  formula,
  data,   
  weights,
  subset, 
  offset, 
  cluster,
  strata,
  na.action = na.pass,
  control = ctree_control(
    teststat = "quad", testtype = "Univ", mincriterion = 0,
    saveinfo = FALSE),
  ytrafo = NULL, 
  scores = NULL, 
  ntree = 500L, 
  perturb = list(replace = FALSE, fraction = 0.632),
  mtry = ceiling(sqrt(nvar)), 
  applyfun = NULL,
  cores = NULL, 
  trace = FALSE,
  weight_variable = NULL) {
  
  ### get the call and the calling environment for .urp_tree
  call <- match.call(expand.dots = FALSE)
  oweights <- NULL
  if (!missing(weights))
    oweights <- weights
  m <- match(c("formula", "data", "subset","weight_variable", "na.action", "offset", "cluster", 
               "scores", "ytrafo", "control"), names(call), 0L)
  ctreecall <- call[c(1L, m)]
  ctreecall$doFit <- FALSE
  if (!is.null(oweights))
    ctreecall$weights <- 1:NROW(oweights)
  ctreecall$control <- control ### put ... into ctree_control()
  ctreecall[[1L]] <- quote(ctree_new)
  tree <- eval(ctreecall, parent.frame())
  
  if (is.null(control$update))
    control$update <- is.function(ytrafo)
  
  d <- tree$d
  updatefun <- tree$update
  d$variable$z <- weight_variable
  nvar <- sum(d$variables$z > 0)
  control$mtry <- mtry
  control$applyfun <- lapply
  
  strata <- d[["(strata)"]]
  if (!is.null(strata)) {
    if (!is.factor(strata)) stop("strata is not a single factor")
  }
  
  probw <- NULL
  iweights <- model.weights(model.frame(d))
  if (!is.null(oweights)) {
    if (is.matrix(oweights)) {
      weights <- oweights[iweights,,drop = FALSE]
    } else {
      weights <- oweights[iweights]
    }
  } else {
    weights <- NULL
  }
  rm(oweights)
  rm(iweights)
  N <- nrow(model.frame(d))
  rw <- NULL
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      if (ncol(weights) == ntree && nrow(weights) == N) {
        rw <- unclass(as.data.frame(weights))
        rw <- lapply(rw, function(w) 
          rep(1:length(w), w))
        weights <- integer(0)
      } else {
        stop(sQuote("weights"), "argument incorrect")
      }
    } else {
      probw <- weights / sum(weights)
    }
  } else {
    weights <- integer(0)
  }
  
  idx <- .start_subset(d)
  if (is.null(rw)) {
    if (is.null(strata)) {
      size <- N
      if (!perturb$replace) size <- floor(size * perturb$fraction)
      rw <- replicate(ntree, 
                      sample(idx, size = size, 
                             replace = perturb$replace, prob = probw),
                      simplify = FALSE)
    } else {
      frac <- if (!perturb$replace) perturb$fraction else 1
      rw <- replicate(ntree, function() 
        do.call("c", tapply(idx, strata, 
                            function(i) 
                              sample(i, size = length(i) * frac, 
                                     replace = perturb$replace, prob = probw[i]))))
    }
  }
  
  ## apply infrastructure for determining split points
  ## use RNGkind("L'Ecuyer-CMRG") to make this reproducible
  if (is.null(applyfun)) {
    applyfun <- if(is.null(cores)) {
      lapply  
    } else {
      function(X, FUN)
        parallel::mclapply(X, FUN, mc.set.seed = TRUE, mc.cores = cores)
    }
  }
  
  trafo <- updatefun(sort(rw[[1]]), integer(0), control, doFit = FALSE)
  if (trace) pb <- txtProgressBar(style = 3) 
  forest <- applyfun(1:ntree, function(b) {
    if (trace) setTxtProgressBar(pb, b/ntree)
    ret <- updatefun(sort(rw[[b]]), integer(0), control)
    # trafo <<- ret$trafo
    ret$nodes
  })
  if (trace) close(pb)
  
  fitted <- data.frame(idx = 1:N)  
  mf <- model.frame(d)
  fitted[[2]] <- mf[, d$variables$y, drop = TRUE]
  names(fitted)[2] <- "(response)"
  if (length(weights) > 0)
    fitted[["(weights)"]] <- weights
  
  ### turn subsets in weights (maybe we can avoid this?)
  rw <- lapply(rw, function(x) as.integer(tabulate(x, nbins = length(idx))))
  
  control$applyfun <- applyfun
  
  ret <- constparties(nodes = forest, data = mf, weights = rw,
                      fitted = fitted, terms = d$terms$all,
                      info = list(call = match.call(), control = control))
  ret$trafo <- trafo
  ret$predictf <- d$terms$z
  class(ret) <- c("cforest", class(ret))
  
  return(ret)
}
ctree_new <- function(formula, data, subset, weights, weight_variable, na.action = na.pass, offset, cluster,
                      control = ctree_control(...), ytrafo = NULL, converged = NULL, scores = NULL,
                      doFit = TRUE, ...) 
{
  
  ## set up model.frame() call
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset","weight_variable", "na.action", "weights",
               "offset", "cluster", "scores"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$yx <- "none"
  if (is.function(ytrafo)) {
    if (all(c("y", "x") %in% names(formals(ytrafo))))
      mf$yx <- "matrix"
  }
  mf$nmax <- control$nmax
  ## evaluate model.frame
  mf[[1L]] <- quote(partykit::extree_data)
  
  d <- eval(mf, parent.frame())
  subset <- .start_subset(d)
  
  weights <- model.weights(model.frame(d))
  
  if (is.function(ytrafo)) {
    if (is.null(control$update))
      control$update <- TRUE
    nf <- names(formals(ytrafo))
    if (all(c("data", "weights", "control") %in% nf))
      ytrafo <- ytrafo(data = d, weights = weights, control = control)
    nf <- names(formals(ytrafo))
    stopifnot(all(c("subset", "weights", "info", "estfun", "object") %in% nf) ||
                all(c("y", "x", "weights", "offset", "start") %in% nf))
  } else {
    if (is.null(control$update))
      control$update <- FALSE
    stopifnot(length(d$variables$x) == 0)
    mfyx <- model.frame(d, yxonly = TRUE)
    mfyx[["(weights)"]] <- mfyx[["(offset)"]] <- NULL
    yvars <- names(mfyx)
    for (yvar in yvars) {
      sc <- d[[yvar, "scores"]]
      if (!is.null(sc))
        attr(mfyx[[yvar]], "scores") <- sc
    }
    Y <- .y2infl(mfyx, response = d$variables$y, ytrafo = ytrafo)
    if (!is.null(iy <- d[["yx", type = "index"]])) {
      Y <- rbind(0, Y)
    }
    ytrafo <- function(subset, weights, info, estfun, object, ...)
      list(estfun = Y, unweighted = TRUE)
    ### unweighted = TRUE prevents estfun / w in extree_fit
  }
  if (is.function(converged)) {
    stopifnot(all(c("data", "weights", "control") %in% names(formals(converged))))
    converged <- converged(d, weights, control = control)
  } else {
    converged <- TRUE
  }
  d$variable$z <- weight_variable
  update <- function(subset, weights, control, doFit = TRUE)
    partykit::extree_fit(data = d, trafo = ytrafo, converged = converged, partyvars = d$variable$z,
                         subset = subset, weights = weights, ctrl = control, doFit = doFit)
  if (!doFit) return(list(d = d, update = update))
  tree <- update(subset = subset, weights = weights, control = control)
  trafo <- tree$trafo
  tree <- tree$nodes
  
  mf <- model.frame(d)
  if (is.null(weights)) weights <- rep(1, nrow(mf))
  
  fitted <- data.frame("(fitted)" = fitted_node(tree, mf),
                       "(weights)" = weights,
                       check.names = FALSE)
  fitted[[3]] <- mf[, d$variables$y, drop = TRUE]
  names(fitted)[3] <- "(response)"
  ret <- party(tree, data = mf, fitted = fitted,
               info = list(call = match.call(), control = control))
  ret$update <- update
  ret$trafo <- trafo
  class(ret) <- c("constparty", class(ret))
  
  ### doesn't work for Surv objects
  # ret$terms <- terms(formula, data = mf)
  ret$terms <- d$terms$all
  ### need to adjust print and plot methods
  ### for multivariate responses
  ### if (length(response) > 1) class(ret) <- "party"
  return(ret)
  
}

rowMatch <- function(A,B) {
  # Rows in A that match the rows in B
  # The row indexes correspond to A
  f <- function(...) paste(..., sep=":")
  if(!is.matrix(B)) B <- matrix(B, 1, length(B))
  a <- do.call("f", as.data.frame(A))
  b <- do.call("f", as.data.frame(B))
  match(b, a)
}

null <- function(n)ifelse(length(n)==0,NA,n)      
addq <- function(x) paste0("`", x, "`")

tree.size <- function( tree) {
  if( is.null( tree)) {
    return( 0)
  } else {
    return( 1 + tree.size( tree$kids[[1]]) + tree.size( tree$kids[[2]]))
  }
}


child_pair_return <- function( x)
{
  if( is.null(x$split))
    c( NA, NA, NA)
  else if( is.null(x$kids[[1]]$split) & is.null(x$kids[[2]]$split))
    c( NA, NA, names(x$info$p.value))
  else if( is.null(x$kids[[1]]$split) & !is.null(x$kids[[2]]$split))
    c( NA,  names(x$kids[[2]]$info$p.value),  names(x$info$p.value))
  else if( !is.null(x$kids[[1]]$split) & is.null(x$kids[[2]]$split))
    c( names(x$kids[[1]]$info$p.value), NA, names(x$info$p.value))
  else
    c( names(x$kids[[1]]$info$p.value), names(x$kids[[2]]$info$p.value), names(x$info$p.value))
}

#var_names <- attr(output$terms, "term.labels")
GetPairVars <- function( t, var_num, var_names)
{  
  #compute the number of nodes in a tree
  nodes_num = tree.size(t)
  
  index_Set <- cbind(1:nodes_num, t$id:(t$id+nodes_num-1))
  #this variable is assigned to check if both the right and left sub tree 
  #from the current node is traversed or not
  # counter[ i] = 0 ===> ith node is not traversed
  # counter[ i] = 1 ===> left subtree of ith node traversed
  # counter[ i] = 2 ===> both left and right subtree traversed.
  counter = rep(0, nodes_num)
  #the following matrix saves the variable at each node and their corresponding left and right child.
  tree_mat = matrix( 0, nrow = nodes_num, ncol = 3)
  colnames(tree_mat)  = c( "left", "right", "parent")
  root1 = root2 = tree_mat2 = list()
  tree_mat3 <- matrix(0, nrow= nodes_num, ncol=var_num+2)
  #tree_mat2 enlists the pair or variables that appear in the tree as child and parent.
  # This list is later utilised to compute the epistasis score, which measures the strength interactions
  # between the variables in the tree
  tree_mat2 = list( list())
  ind = 1
  root1[[ ind]] = node = t
  #nodeID indicates the index of the node under study. Note that since neither left nor right subtree is encountered the counter is set 0
  counter[match(node$id,index_Set[,2])] = 0
  #returns the children of each node in the corresponding row of the matrix.
  tree_mat[match(node$id,index_Set[,2]), ] = child_pair_return( node)
  checkNull <- function( x)ifelse( is.null(x$split), NA, names(x$info$p.value))
  VarList <- function( x, y) c( checkNull( x), checkNull( y))
  colnames(tree_mat3)=c("Immediate_parent",var_names, "Current_node")
  tree_mat3[node$id, 1] = "root"
  tree_mat3[node$id,var_num+2] = checkNull(node)
  tree_mat3[node$id, match(colnames(node$info$criterion), var_names)+1]=1
  #tree_mat3[node$id,var_num+2] = paste("V",node$split$varid-1, sep="")
  #tree_mat2[[ paste( "V", node$split$varid-1,sep="")]] = rbind( tree_mat2[[ paste( "V", node$split$varid-1,sep="")]], VarList( node, root1[[ ind]]))
  # The tree_mat matrix reports NA if a node is a terminal node, i.e. tree_mat[ i, 2] = NA means left child of node i is empty and
  # tree_mat[ i, 3] = NA means right child of node i is empty
  # Hence if the tree_mat have any entry 0 (not assigned) means the tree is not fully traversed.
  #tree_root <- node
  
  
  while( any( na.omit( c( tree_mat))==0))
  {
    # to traverse the left subtree
    # !( node$terminal) ~~ The node cannot be a terminal node and 
    # counter[ node$nodeID]!=2 ~~ both the left subtree and right subtree of the node cannot be traveresed 
    while( !( is.null(node$split)) & counter[ match(node$id, index_Set[,2])]!=2)
    {
      tree_mat3[node$kids[[1]]$id, 1] =  checkNull(node)
      node = node$kids[[1]]
      tree_mat3[node$id, match(colnames(node$info$criterion),var_names)+1]=1
      tree_mat3[node$id,var_num+2] =  checkNull(node)
      counter[ match(root1[[ ind]]$id, index_Set[,2])] = counter[ match(root1[[ ind]]$id, index_Set[,2])] + 1
      tree_mat[ match(node$id, index_Set[,2]), ] = child_pair_return( node)
      #if(!(is.na(tree_mat[ match(node$id, index_Set[,2]), 3])))
        #if(tree_mat[ match(node$id, index_Set[,2]), 3]=="V")
          #tree_mat[ match(node$id, index_Set[,2]), 3]=NA
      ind2 = ind
      if( ind==1)
      {
        if( !any( is.na( VarList( node, root1[[ ind]]))))
          #if(root1[[ ind]]$id!=1)
          tree_mat2[[ names(node$info$p.value)]] = rbind( tree_mat2[[ names(node$info$p.value)]], VarList( node, root1[[ ind]]))
      }
      else
      {
        while( ind2!=0)
        { root2 = root1[[ ind2]]
        if( !any( is.na( VarList( node, root2))))
          #if(root2$id!=1)
          tree_mat2[[ names(node$info$p.value)]] = rbind( tree_mat2[[ names(node$info$p.value)]], VarList(node, root2))
        ind2 = ind2-1
        }
      }
      
      ind = ind+1
      root1[[ ind]] = node
      
    }
    # to traverse the right subtree
    ind=ind-1
    #tree_mat3[node$id, 1] = paste( "V",  node$split$varid-1, sep="")
    node=root1[[ ind]]
    #tree_mat3[node$id, match(colnames(node$info$criterion), paste("SNP", 1:var_num, sep="_"))+1]=1
    #tree_mat3[node$id,1002] = paste("V",node$split$varid-1, sep="")
    if(counter[match( node$id, index_Set[,2])]!=2)
    {
      counter[ match( node$id, index_Set[,2])] = counter[match( node$id, index_Set[,2])] + 1
      tree_mat3[node$kids[[2]]$id, 1] =  checkNull(node)
      node = node$kids[[2]]
      tree_mat3[node$id, match(colnames(node$info$criterion), var_names)+1]=1
      tree_mat3[node$id,var_num+2] =  checkNull(node)
      tree_mat[ match( node$id, index_Set[,2]), ] = child_pair_return( node)
      #if(!(is.na(tree_mat[ match(node$id, index_Set[,2]), 3])))
       # if(tree_mat[ match(node$id, index_Set[,2]), 3]=="V")
        #  tree_mat[ match(node$id, index_Set[,2]), 3]=NA
      
      ind2 = ind
      if( ind==1)
      {
        if( !any( is.na( VarList( node, root1[[ ind]]))))
          #if(root1[[ ind]]$id!=1)
          tree_mat2[[ names(node$info$p.value)]] = rbind( tree_mat2[[ names(node$info$p.value)]], VarList( node, root1[[ ind]]))
      }
      else
      {
        while( ind2!=0)
        { root2 = root1[[ ind2]]
        if( !any( is.na( VarList( node, root2))))
          #if(root2$id!=1)
          tree_mat2[[ names(node$info$p.value)]] = rbind( tree_mat2[[ names(node$info$p.value)]], VarList(node, root2))
        ind2 = ind2-1
        }
      }
      ind = ind + 1
      root1[[ ind]] = node
      
      
      
    }
  }
  #tree_mat2[[tree_root$psplit$variableName]] <- cbind(tree_root$psplit$variableName, setdiff(unique(na.omit(c(tree_mat))),tree_root$psplit$variableName))
  return(list(tree_mat3, tree_mat2, tree_mat))
}

rowmatch <- function( A, B) { 
  # Rows in A that match the rows in B
  f <- function(...) paste(..., sep=":")
  if( !is.matrix(B)) B <- matrix( B, 1, length( B))
  a <- do.call( "f", as.data.frame( A))
  b <- do.call( "f", as.data.frame( B))
  ifelse( is.na( match( b, a)), 0, match( b, a))
}

row.matches <- function(y, X){
  i <- seq(nrow(X))
  j <- 0
  while(length(i) && (j <- j + 1) <= ncol(X))
    i <- i[X[i, j] == y[j]]
  i
}

swap <-function(vec)
{
  temp <- vec[1]
  vec[1] <- vec[2]
  vec[2] <- temp
  return(vec)  
}


Var_Interaction_Proportions <- function( tree_mat3)
{
  #var_num=ncol(tree_mat3)-2
  #library(varhandle)
  #tree_mat3 <- as.data.frame(tree_mat3) 
  #colnames(tree_mat3)[2:(var_num+1)] = c(paste("SNP", 1:var_num, sep=""))
  #tree_mat3[,1:(ncol(tree_mat3)-1)] <- sapply(1:(ncol(tree_mat3)-1), function(i)unfactor(tree_mat3[,i]))
  
  #tree_mat3[is.na(tree_mat3)]=""
  #pos_to_check = which(tree_mat3[,"Current_node"]=="")
  #pos_to_check = pos_to_check[which(tree_mat3[pos_to_check-1, "Current_node"]!=""))]
  #pos_to_check = na.omit(sapply(pos_to_check, function(pos)ifelse(all(tree_mat3[pos,c("Current_node", "Immediate_parent")]==tree_mat3[pos+1,c("Current_node", "Immediate_parent")]), pos,NA)))
  #pos_to_check = setdiff(pos_to_check,2)
  
  var_num=ncol(tree_mat3)-2
  library(varhandle)
  tree_mat3 <- as.data.frame(tree_mat3) 
  tree_mat3[,1:(ncol(tree_mat3)-1)] <- sapply(1:(ncol(tree_mat3)-1), function(i)as.matrix(tree_mat3[,i]))
  
  
  pos_to_check = which(is.na(tree_mat3[,"Current_node"]))
  pos_to_check = pos_to_check[which(!is.na(tree_mat3[pos_to_check-1, "Current_node"]))]
  pos_to_check = setdiff(pos_to_check,2)
  
  pos_start=1
  interaction <- list()
  ind=1
  if(length(pos_to_check)!=0)
  {
    for(pos in pos_to_check)
    {
      root_of_age <- tree_mat3[ pos, "Immediate_parent"]
      interaction[[ind]] <- vector()
      pos_end = pos
      
      while(root_of_age!="root"&length(pos_end)!=0)
      {
        interaction[[ind]] = append(interaction[[ind]], root_of_age)
        pos_new = which(tree_mat3[pos_start:(pos_end-1), "Current_node"]==tree_mat3[pos_end, "Immediate_parent"])
        root_of_age = tree_mat3[tail(pos_new,1),"Immediate_parent"]
        pos_end=tail(pos_new,1)
      }
      ind=ind+1
      #l=l+1
    }
  }
  
  return(interaction)
  
}

getInteractionMatrix <- function(output, Snp_List)
{
  
  ntree=length(output$nodes)
  var_num=ncol(output$data)-1
  var_names <- attr(output$terms, "term.labels" )
  Interaction_Score1 <- Interaction_Score2 <- var_interaction_score <- var_interaction_score_buf<- matrix(0, nrow=var_num, ncol=var_num)
  colnames(Interaction_Score1) <- rownames(Interaction_Score1) <- var_names
  colnames(Interaction_Score2) <- rownames(Interaction_Score2) <- var_names
  
   
  Interaction_Lists <- list()
  #Snp_List <- c(SNPs_Selected, SNPs_for_cf)
  #Snp_Map <- cbind(Snp_List, paste("V", 1:length(Snp_List), sep=""))
  #paste("SNP_", 1:10, sep="")  
  Var_Interaction_Mat <- list()
  for( ind in 1:ntree)
  {
    tree_ind =output$nodes[[ind]]
    #print(length(tree_ind$kids))
    Interaction_Lists[[ind]] <-list()
    
    if(length(tree_ind$kids)==0)
      Var_Interaction_Mat[[ind]]=list()
    else
    { 
      Var_Interaction <- GetPairVars(tree_ind, var_num, var_names )
      Var_Proportion <- Var_Interaction_Proportions(Var_Interaction[[1]])
       if(length(Var_Proportion)==0)
         Interaction_Lists[[ind]] <- list()
       else
         Interaction_Lists[[ind]] <- Var_Proportion
      if(length(unlist(Var_Interaction[[2]]))!=0)
      {
         Var_Interaction_Mat[[ind]] <- matrix( unlist( do.call( "rbind", Var_Interaction[[2]] )), ncol = 2)
        # replacing the SNP variable by the original SNP name
         Var_Interaction_Mat[[ind]] <-  t(apply(Var_Interaction_Mat[[ind]], 1, grr::sort2))
      }
    else 
      Var_Interaction_Mat[[ind]]=list()
    }
    for(i_l in 1:length(Interaction_Lists))
      Interaction_Lists[[i_l]] <- sapply(Interaction_Lists[[i_l]], function(list_member)unique(list_member))
  }
  

  
  len_mat <- sapply(1:length(Var_Interaction_Mat),function(i)length(Var_Interaction_Mat[[i]]))
  
  if(any(len_mat!=0))
    # This lists the unique interaction in the entire forest
  {
    Interaction <- matrix( unlist( do.call( "rbind", Var_Interaction_Mat)), ncol=2)
    Unique_Interaction <- unique( Interaction)
  
  
  # Importance of each SNP interacting with age, this score in a way measures the strength of an SNP is the
  # aging of the phenotype
    Importance_Score <- cbind( Variable1 = Unique_Interaction[ ,1], Variable2 = Unique_Interaction[ ,2], Score = rep( 0, nrow(Unique_Interaction)))
    Importance_Score[,1:2] <- t(apply(Importance_Score[,1:2], 1, grr::sort2))
    Importance_Score <- unique(Importance_Score)
    system.time(for(ind in 1:nrow( Importance_Score))
    Importance_Score[ind, 3]<- sum(sapply(1:ntree, function(j)ifelse(length(Var_Interaction_Mat[[j]])!=0,
                                                                     {ifelse(is.na(row.match( Importance_Score[ind, 1:2], Var_Interaction_Mat[[j]])),0,1)}, 0))/ntree))
    if(nrow(Importance_Score)>1)
    {
      Importance_Score <- as.data.frame( Importance_Score[ order( Importance_Score[, 3], decreasing = TRUE), ])
      
      Importance_Score <- Importance_Score %>% filter(Score>1/ntree)
    #g<- graph_from_edgelist( as.matrix(Importance_Score[,1:2]), directed=F)
    }
  }
  else
    Importance_Score =NA
  
  for(ind in 1:length(Interaction_Lists))
  {
    n_I=length(Interaction_Lists[[ind]])
    if(n_I>0)
      Interaction_Lists[[ind]] <- lapply(1:n_I, function(k)interaction_list(Interaction_Lists[[ind]][[k]]))
    #Interaction_Lists[[ind]] <- lapply(1:n_I,function(k){lapply(3:(length(Interaction_Lists[[ind]][[k]])), function(j)Interaction_Lists[[ind]][[k]][1:j])})
    
  }
  All_Interactions <- lapply(1:length(Interaction_Lists), 
                             function(tree)
                             {if(length(Interaction_Lists[[tree]])!=0)
                             {lapply(1:length(Interaction_Lists[[tree]]), function(k)
                             {if(length(Interaction_Lists[[tree]][[k]])>0)
                               do.call(qpcR:::rbind.na,Interaction_Lists[[tree]][[k]])
                               else return(NA)})}})
  All_Interactions_Full <- lapply(1:length(All_Interactions), function(k)if(!is.null(All_Interactions[[k]])){if(length(All_Interactions[[k]])==1)return(All_Interactions[[k]][[1]])else return(do.call(qpcR:::rbind.na, All_Interactions[[k]]))})
  
  All_Interactions <- do.call(qpcR:::rbind.na,  All_Interactions_Full)
  if(ncol(All_Interactions)==0)
    return(list( Interaction_Score1, All_Interactions))
  else
  {
    colnames(All_Interactions) <- paste("Node", 1:ncol(All_Interactions), sep="")
    All_Interactions <- unique(All_Interactions)
    #All_Interactions <- All_Interactions%>% distinct()
    row_all_na <- sapply(1:nrow(All_Interactions), function(row)ifelse( all(is.na(All_Interactions[row,])), row,NA))
    if(length(na.omit(row_all_na))>0)
      All_Interactions <-All_Interactions[-na.omit(row_all_na),]
    
    #All_Interactions <- All_Interactions[-which(is.na(All_Interactions[,2])),]
    
    All_Interactions[is.na(All_Interactions)]="z"
    All_Interactions[All_Interactions==""]="z"
    All_Interactions_New <- t(apply(All_Interactions, 1, grr::sort2))
    All_Interactions_New <-  unique(All_Interactions_New)
    All_Interactions_New[All_Interactions_New=="z"]=""
    All_Interactions_Stats_New <- cbind( All_Interactions_New, matrix(0, ncol=ntree, nrow=nrow(All_Interactions_New), dimnames = list(paste("Interaction_",1:nrow(All_Interactions_New), sep=""),c( paste("Tree_", 1:ntree, sep="")))))
    
    
    system.time(for(i in  1:length(Interaction_Lists))
    {
      if(!is.null(nrow(All_Interactions_Full[[i]]))&!all(is.na(All_Interactions_Full[[i]])))
      {
        All_Interactions_Full[[i]] <- as.matrix(All_Interactions_Full[[i]])
        colnames(All_Interactions_Full[[i]])<-paste("V",1:ncol(All_Interactions_Full[[i]]), sep="_")
        All_Interactions_Full[[i]][is.na(All_Interactions_Full[[i]])]="z"
        All_Interactions_Full[[i]][All_Interactions_Full[[i]]==""]="z"
        
        All_Interactions_Full[[i]] <-  t(apply(All_Interactions_Full[[i]], 1, grr::sort2))
        All_Interactions_Full[[i]][All_Interactions_Full[[i]]=="z"]=""
        #occ <- sapply(1:nrow(All_Interactions_New), function(interaction) sum(!is.na(row.match(All_Interactions_New[interaction,1:ncol(All_Interactions_Full[[i]])], All_Interactions_Full[[i]]))))
        occ <-  rowMatch(All_Interactions_Full[[i]], All_Interactions_New[,1:ncol(All_Interactions_Full[[i]])])
        occ[!is.na(occ)]=1
        occ[is.na(occ)]=0
        All_Interactions_Stats_New[,paste("Tree_", i, sep="")] =occ
      }
    })
    
    All_Interactions_Stats_New <- as.data.frame(All_Interactions_Stats_New)
    All_Interactions_Stats_New[, paste( "Tree_",1:ntree, sep="")]<- lapply( All_Interactions_Stats_New[, paste( "Tree_",1:ntree, sep="")], as.numeric)
    All_Interactions_Stats_New$Forest_Score <- rowSums( All_Interactions_Stats_New[, paste( "Tree_",1:ntree, sep="")])
    All_Interactions_Stats_New <- All_Interactions_Stats_New[ order( All_Interactions_Stats_New$Forest_Score),]
    All_Interactions_Stats_New <- All_Interactions_Stats_New[,-match(paste("Tree_",1:ntree,sep=""), colnames(All_Interactions_Stats_New))]
    All_Interactions_Stats_New <- All_Interactions_Stats_New[- which(All_Interactions_Stats_New[,2]==""),]
    
    
    return( list( Importance_Score, All_Interactions_Stats_New))
  }
}

GenerateInteractionList <- function(Interaction_List, Importance_Score_list)
{
  
  All_Interactions_Stats <- do.call(qpcR:::rbind.na, Interaction_List)
  All_Interactions_Stats <-All_Interactions_Stats %>% distinct()
  All_Interactions_Stats[All_Interactions_Stats==""]="z"
  All_Interactions_Stats[is.na(All_Interactions_Stats)]="z"
  All_Interactions_Stats <- t(apply(All_Interactions_Stats, 1, grr::sort2))
  All_Interactions_Stats[All_Interactions_Stats=="z"]=""
  colnames(All_Interactions_Stats) <- paste("Node",1:ncol(All_Interactions_Stats), sep="")
  All_Interactions_Stats <- as.data.frame(All_Interactions_Stats)
  All_Interactions_Stats <-All_Interactions_Stats %>% distinct()
  NiterChild=10

  All_Interactions <- All_Interactions_Stats
  All_Interactions_Stats <- cbind(All_Interactions_Stats, matrix(0, ncol=NiterChild, nrow=nrow(All_Interactions_Stats), dimnames = list(paste("Interaction_",1:nrow(All_Interactions_Stats), sep=""),c( paste("Forest_", 1:NiterChild, sep="")))))
  
  system.time(for(i in 1:NiterChild)
  {  
    mat_to_check <- Importance_Score_list[[i]][[2]]
    mat_to_check <- mat_to_check %>% filter(Forest_Score>0)
    if(length(which(colSums(mat_to_check=="")==nrow(mat_to_check)))>0)
      mat_to_check <- mat_to_check[, which(colSums(mat_to_check=="")==nrow(mat_to_check))]
    
    system.time(pos_match <- matchRow(nrow(mat_to_check), as.matrix(All_Interactions_Stats[,1:(ncol(mat_to_check)-1)]), as.matrix(mat_to_check[, -ncol(mat_to_check)])))
    #system.time(pos_match <-vapply(1:nrow(mat_to_check), function(row)which(compareToRow(as.matrix(All_Interactions_Stats[,1:(ncol(mat_to_check)-1)]),as.vector(unlist(mat_to_check[row,-ncol(mat_to_check)])))), numeric(1)))
    
    All_Interactions_Stats[pos_match,paste("Forest_",i,sep="")] = mat_to_check$Forest_Score
    #All_Interactions_Stats[[paste("Forest_",i,sep="")]][All_Interactions_Stats[[paste("Forest_",i,sep="")]]>1]=1
  })
  
  All_Interactions_Stats$Sum_Forest_Score <- sapply( 1:nrow(All_Interactions_Stats), function(i) sum( as.numeric(All_Interactions_Stats[ i, paste("Forest_", 1:10,sep="")]), na.rm=TRUE) )
  All_Interactions_Stats$Median_Forest_Score <- rowMedians(as.matrix(All_Interactions_Stats[,paste("Forest",1:10, sep="_")]))
  
  All_Interactions_Stats<- All_Interactions_Stats%>% distinct()
  All_Interactions_Stats<- All_Interactions_Stats[order(All_Interactions_Stats$Sum_Forest_Score),]
  #All_Interactions_Stats$Rank <- rank(max(All_Interactions_Stats$Median_Forest_Score)-All_Interactions_Stats$Median_Forest_Score)
  
  return(All_Interactions_Stats)
}

plotSNPInteraction <- function(data_epistasis, snp_list)
{
  
  data_epistasis$Genotype_Combination <- factor(sapply(1:nrow(data_epistasis),FUN=function(x)paste(unfactor(data_epistasis[x,snp_list]),collapse = "")))
  g <- list()
  xlabs <- paste(levels(data_epistasis$Genotype_Combination),"\n (N=",table(data_epistasis$Genotype_Combination),")",sep="")
  g[[1]] <- ggplot(data_epistasis, aes(x=Genotype_Combination, y=PHENOTYPE, fill=Genotype_Combination))+scale_fill_viridis(discrete = TRUE, alpha=0.6)+scale_x_discrete(labels=xlabs)+
    geom_jitter(color="black", size=0.4, alpha=0.9) + geom_boxplot()+theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size=11))
  for(i in 1:length(snp_list))
  {
   xlabs <- paste(levels(data_epistasis[[snp_list[i]]]),"\n (N=",table(data_epistasis[[snp_list[i]]]),")",sep="")
   g[[i+1]] <- ggplot(data_epistasis, aes(x=get(snp_list[i]), y=PHENOTYPE, fill=get(snp_list[i]))) +scale_fill_viridis( discrete = TRUE, alpha=0.1)+scale_x_discrete(labels=xlabs) +
    geom_boxplot()+theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size=11))+xlab(snp_list[i])
  }
  fit1 <- lm(as.formula(paste("PHENOTYPE~",paste(addq(snp_list), collapse="+"))), data= data_epistasis)
  fit2 <- lm(as.formula(paste("PHENOTYPE~Genotype_Combination")), data= data_epistasis)
  an <- anova(fit1, fit2, test="Chisq")
  lay <- rbind(2:(length(snp_list)+1),rep(1, length(snp_list)))
  #g[[i+2]]<- grid.table(an)%>%kable()
  grid.arrange(grobs=g, layout_matrix = lay,  top = textGrob(paste("Testing the effect of the n-way interaction between", paste(snp_list, collapse = ","), sep=""),gp=gpar(fontsize=10,font=2)))

  print(an)
 }


RowContain <- function(set, mat) sapply(1:nrow(mat), function(r)ifelse(all(set %in% mat[r,]), r, NA))

interaction_list <- function(var)
{
  if(length(var)>1)
  {
    bin_mat <- eval(parse(text=paste("expand.grid(", paste(rep("0:1",length(var)), collapse=","),")", sep="")))
    bin_mat <- bin_mat[-1,]
    interactions <- lapply(1:nrow(bin_mat), function(row)c(var[which(bin_mat[row,]==1)]))
    return(interactions)
  }
  else
    return(list())
}

library(Rcpp)
cppFunction(
'NumericVector matchRow(int I, CharacterMatrix x, CharacterMatrix y) {
   int nrow = x.nrow();
  NumericVector out(nrow);
  LogicalVector s(nrow);

  for (int i = 0; i < I; i++) {
   s =compareToRow(x, y.row(i));
   for(int ind = 0; ind < s.size(); ind++)
    {
        if (s[ind] == true)
        {
            out[i] = ind+1;
            break;
        }
    }
  
  }
  return out;
}', includes= 'LogicalVector compareToRow(CharacterMatrix x, CharacterVector y) {
  const int nr = x.nrow();
  const int nc = x.ncol();
  LogicalVector ret(nr, true);
  for (int j=0; j < nr; ++j) {
    for (int k=0; k < nc; ++k) {
      if (x(j, k) != y[k]) {
        ret[j] = false;
        break;
      }
    }
  }
  return ret;
}')




#library(Rcpp)
#cppFunction(
#  'LogicalVector compareToRow(CharacterMatrix x, CharacterVector y) {
#  const int nr = x.nrow();
#  const int nc = x.ncol();
#  LogicalVector ret(nr, true);
#  for (int j=0; j < nr; ++j) {
#    for (int k=0; k < nc; ++k) {
#      if (x(j, k) != y[k]) {
#        ret[j] = false;
#        break;
#      }
#    }
#  }
#  return ret;
#}')

#Test Run
#system.time(pos_match <-vapply(1:nrow(mat_to_check), function(row)which(compareToRow(as.matrix(All_Interactions_Stats[,1:(ncol(mat_to_check)-1)]),as.vector(unlist(mat_to_check[row,-ncol(mat_to_check)])))), numeric(1)))
#Takes 236 secs
#system.time(pos_match2 <- matchRow(10, as.matrix(All_Interactions_Stats[,1:(ncol(mat_to_check)-1)]), as.matrix(mat_to_check[, -ncol(mat_to_check)])))
#Takes 57 secs




Max.Test.Age1 <- function(Cluster)
{
  
  #For each combination we need to find the beta that is the effect due to the age or the temporal effect
  # Y~\betaxt + \alpha + \epsilon
  SNP_Level_Summary <- table(Cluster$SNP_Clubbed)
  #if(length(which(SNP_Level_Summary<10))>0)
   #Cluster <- Cluster%>% filter(!(SNP_Clubbed%in%names(which(SNP_Level_Summary<10))))
  
  Cluster$SNP_Clubbed <- factor(Cluster$SNP_Clubbed)
  SNP_levels <- levels(Cluster$SNP_Clubbed)
  reference <- SNP_levels[1]
  Beta_Estimate <-Beta_SE <- rep(0, length(SNP_levels))
  T <- Stderror<-  rep(0, (nlevels(Cluster$SNP_Clubbed)-1))
  df <- rep(0, length(SNP_levels)-1)
  
  for(snp_level in setdiff(SNP_levels, reference))
  {
  j=match(snp_level,SNP_levels)-1
   Cluster_Of_Interest <- Cluster %>% filter(SNP_Clubbed%in% c(reference, snp_level))
    fit1 <- t.test(value~SNP_Clubbed, data=Cluster_Of_Interest)
    T[j]=fit1$statistic
    Beta_Estimate[j+1]=fit1$estimate[2]
    Beta_Estimate[1]=fit1$estimate[1]
    df[j]=fit1$parameter
    Stderror[j]=fit1$stderr
  }
  
  Ti.Table <- data.frame(diff_combi=Beta_Estimate[-1]-Beta_Estimate[1], Stderror= Stderror, estimate=T, DF=df )
  rownames( Ti.Table) <- setdiff(SNP_levels, reference)
  Ti.Table <- Ti.Table%>%mutate(pval.unadj= 1-(pt(abs(T), DF)-pt(-abs(T), DF)))
  Ti.Table <- Ti.Table%>%mutate(pval.adj= 1-(pt(max(abs(T)), DF)-pt(-max(abs(T)), DF)))
  R=diag(nlevels(Cluster$SNP_Clubbed)-1) 
  R[upper.tri(R) | lower.tri(R)] = 0.5
  pval <-1- pmvnorm(lower=rep(-max( abs(Ti.Table$estimate)),length(SNP_levels)-1), upper=rep(max(abs(Ti.Table$estimate)),length(SNP_levels)-1), mean=rep(0, length(SNP_levels)-1)
                     , sigma=R, algorithm = GenzBretz())
  return( pval)
}
