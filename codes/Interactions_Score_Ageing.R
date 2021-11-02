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
#### Packages/libraries ####

if(!require(forcats)){
  install.packages("forcats")
  library(forcats)
}

if(!require(dplyr)){
  library(partykit)
}

if(!require(ggplot2)){
  library(ggplot2)
}

if(!require(stringr)){
  library(stringr)
}

if(!require(gridExtra)){
  install.packages("gridExtra")
  library(gridExtra)
}

if(!require(dplyr)){
  library(dplyr)
}

# Plotting random clusters 
library(gridExtra)
library(grid)
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

null <- function(n)ifelse(length(n)==0,NA,n)

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
                      doFit = TRUE, ...) {
  
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
    extree_fit(data = d, trafo = ytrafo, converged = converged, partyvars = d$variable$z,
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

get_cTree <- function(cf, k=1) {
  dt <- cf@data@get("input")
  tr <- party:::prettytree(cf@ensemble[[k]], names(dt))
  tr_updated <- update_tree(tr, dt)
  new("BinaryTree", tree=tr_updated, data=cf@data, responses=cf@responses, 
      cond_distr_response=cf@cond_distr_response, predict_response=cf@predict_response)
}

update_tree <- function(x, dt) {
  x <- update_weights(x, dt)
  if(!x$terminal) {
    x$left <- update_tree(x$left, dt)
    x$right <- update_tree(x$right, dt)   
  } 
  x
}

update_weights <- function(x, dt) {
  splt <- x$psplit
  spltClass <- attr(splt,"class")
  spltVarName <- splt$variableName
  spltVar <- dt[,spltVarName]
  spltVarLev <- levels(spltVar)
  if (!is.null(spltClass)) {
    if (spltClass=="nominalSplit") {
      attr(x$psplit$splitpoint,"levels") <- spltVarLev   
      filt <- spltVar %in% spltVarLev[as.logical(x$psplit$splitpoint)] 
    } else {
      filt <- (spltVar <= splt$splitpoint)
    }
    x$left$weights <- as.numeric(filt)
    x$right$weights <- as.numeric(!filt)
  }
  x
}

tree.size <- function( tree) {
  if( is.null( tree)) {
    return( 0)
  } else {
    return( 1 + tree.size( tree$kids[[1]]) + tree.size( tree$kids[[2]]))
  }
}


child_pair_return <- function( x, var_names)
{
  #var_names <- c("age", var_names)
  if( is.null(x$kids))
    c( NA, NA, NA)
  else if( is.null(x$kids[[1]]$split) & is.null(x$kids[[2]]$split))
    c( NA, NA, var_names[x$split$varid-1])
  else if( is.null(x$kids[[1]]$split) & !is.null(x$kids[[2]]$split))
    c( NA, var_names[x$kids[[2]]$split$varid-1], var_names[x$split$varid-1])
  else if( !is.null(x$kids[[1]]$split) & is.null(x$kids[[2]]$split))
    c( var_names[x$kids[[1]]$split$varid-1], NA,  var_names[x$split$varid-1])
  else
    c( var_names[x$kids[[1]]$split$varid-1], var_names[x$kids[[2]]$split$varid-1], var_names[x$split$varid-1])
}

GetPairVars <- function( t, var_num, var_names)
{  
  #compute the number of nodes in a tree
  nodes_num = tree.size(t)
  node_summary = nodeapply(t, ids = nodeids(t), FUN = function(n) n$split$varid)
  var_names_selected = sapply(1:nodes_num, function(k) {if(is.null(node_summary[[k]]))node_summary[[k]]else{ifelse(node_summary[[k]]==2, "age", var_names[node_summary[[k]]-2])}})
  
  index_Set <- cbind(1:nodes_num, var_names_selected)
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
  tree_mat3 <- matrix(0, nrow= nodes_num, ncol=var_num+3)
  #tree_mat2 enlists the pair or variables that appear in the tree as child and parent.
  # This list is later utilised to compute the epistasis score, which measures the strength interactions
  # between the variables in the tree
  tree_mat2 = list( list())
  ind = 1
  root1[[ ind]] = node = t
  #nodeID indicates the index of the node under study. Note that since neither left nor right subtree is encountered the counter is set 0
  counter[match(node$id,index_Set[,1])] = 0
  #returns the children of each node in the corresponding row of the matrix.
  tree_mat[match(node$id,index_Set[,1]), ] = child_pair_return( node, var_names)
  var_names <- c("age", var_names)
  checkNull <- function( x)ifelse( is.null(x$split), NA, var_names[x$split$varid-1])
  VarList <- function( x, y) c( checkNull( x), checkNull( y))
  colnames(tree_mat3)=c("Immediate_parent",var_names,"Current_node")
  tree_mat3[node$id, 1] = "root"
  tree_mat3[node$id,var_num+3] = checkNull(node)
  tree_mat3[node$id, match(colnames(node$info$criterion), var_names)+1]=1
  #tree_mat3[node$id,var_num+2] = paste("SNP_",node$split$varid-1, sep="")
  #tree_mat2[[ paste( "SNP_", node$split$varid-1,sep="")]] = rbind( tree_mat2[[ paste( "SNP_", node$split$varid-1,sep="")]], VarList( node, root1[[ ind]]))
  # The tree_mat matrix reports NA if a node is a terminal node, i.e. tree_mat[ i, 2] = NA means left child of node i is empty and
  # tree_mat[ i, 3] = NA means right child of node i is empty
  # Hence if the tree_mat have any entry 0 (not assigned) means the tree is not fully traversed.
  #tree_root <- node
  
  while( any( na.omit( c( tree_mat))==0))
  {
    # to traverse the left subtree
    # !( node$terminal) ~~ The node cannot be a terminal node and 
    # counter[ node$nodeID]!=2 ~~ both the left subtree and right subtree of the node cannot be traveresed 
    while( !( is.null(node$split)) & counter[ match(node$id, index_Set[,1])]!=2)
    {
      
      tree_mat3[node$kids[[1]]$id, 1] =  checkNull(node)
      node = node$kids[[1]]
      tree_mat3[node$id, na.omit(match(colnames(node$info$criterion), var_names))+1]=1
      tree_mat3[node$id,var_num+3] =  checkNull(node)
      counter[ match(root1[[ ind]]$id, index_Set[,1])] = counter[ match(root1[[ ind]]$id, index_Set[,1])] + 1
      tree_mat[ match(node$id, index_Set[,1]), ] = child_pair_return( node, var_names)
      if(!(is.na(tree_mat[ match(node$id, index_Set[,1]), 3])))
        if(tree_mat[ match(node$id, index_Set[,1]), 3]=="SNP_")
          tree_mat[ match(node$id, index_Set[,1]), 3]=NA
      ind2 = ind
      if( ind==1)
      {
        if( !any( is.na( VarList( node, root1[[ ind]]))))
          #if(root1[[ ind]]$id!=1)
          tree_mat2[[var_names[node$split$varid-1]]] = rbind( tree_mat2[[ var_names[node$split$varid-1]]], VarList( node, root1[[ ind]]))
      }
      else
      {
        while( ind2!=0)
        { root2 = root1[[ ind2]]
        if( !any( is.na( VarList( node, root2))))
          #if(root2$id!=1)
          tree_mat2[[var_names[node$split$varid-1]]] = rbind( tree_mat2[[var_names[node$split$varid-1]]], VarList(node, root2))
        ind2 = ind2-1
        }
      }
      
      ind = ind+1
      root1[[ ind]] = node
      
    }
    # to traverse the right subtree
    ind=ind-1
    #tree_mat3[node$id, 1] = paste( "SNP_",  node$split$varid-1, sep="")
    node=root1[[ ind]]
    #tree_mat3[node$id, match(colnames(node$info$criterion), paste("SNP", 1:var_num, sep="_"))+1]=1
    #tree_mat3[node$id,1002] = paste("SNP_",node$split$varid-1, sep="")
    if(counter[match( node$id, index_Set[,1])]!=2)
    {
      counter[ match( node$id, index_Set[,1])] = counter[match( node$id, index_Set[,1])] + 1
      tree_mat3[node$kids[[2]]$id, 1] =  checkNull(node)
      node = node$kids[[2]]
      tree_mat3[node$id, na.omit(match(colnames(node$info$criterion), var_names))+1]=1
      tree_mat3[node$id,var_num+3] =  checkNull(node)
      tree_mat[ match( node$id, index_Set[,1]), ] = child_pair_return( node, var_names)
      if(!(is.na(tree_mat[ match(node$id, index_Set[,1]), 3])))
        if(tree_mat[ match(node$id, index_Set[,1]), 3]=="SNP_")
          tree_mat[ match(node$id, index_Set[,1]), 3]=NA
      
      ind2 = ind
      if( ind==1)
      { if( !any( is.na( VarList( node, root1[[ind]]))))
        #if(root1[[ ind]]$id!=1)
        tree_mat2[[ var_names[node$split$varid-1]]] = rbind( tree_mat2[[ var_names[node$split$varid-1]]], VarList( node, root1[[ ind]]))
      }
      else
      {
        while( ind2!=0)
        { root2 = root1[[ ind2]]
        if(!any( is.na( VarList( node, root2))))
          #if(root2$id!=1)
          tree_mat2[[ var_names[node$split$varid-1]]] = rbind( tree_mat2[[  var_names[node$split$varid-1]]], VarList( node, root2))
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
  var_num=ncol(tree_mat3)-2
  library(varhandle)
  tree_mat3 <- as.data.frame(tree_mat3) 
  tree_mat3[,1:(ncol(tree_mat3)-1)] <- sapply(1:(ncol(tree_mat3)-1), function(i)as.matrix(tree_mat3[,i]))
  var_sampling_matrix <- var_sampling_matrix_true<-matrix(0, nrow=0, ncol=var_num+2)
  colnames(var_sampling_matrix)=colnames(tree_mat3)
  
  interaction <- list()
  count=1
  pos_start=1
  positions=which(tree_mat3$Current_node=="age")
  vsm_row=0
  rowSumBinary <- function(a,b)
  {
    c=as.numeric(a)+as.numeric(b)
    c[c>1]=1
    return (c)
  }
  for(pos in positions)
  {
    root_of_age = tree_mat3$Immediate_parent[pos]
    interaction[[count]]="age"
    pos_end=pos
    ind=1
    interaction_summary <- matrix(0, nrow=0, ncol=var_num)
    while(root_of_age!="root")
    {
      interaction[[count]] = append(interaction[[count]], root_of_age)
      pos_new = which(tree_mat3$Current_node[pos_start:(pos_end-1)]==tree_mat3$Immediate_parent[pos_end])
      root_of_age = tree_mat3$Immediate_parent[tail(pos_new,1)]
      #if(root_of_age=="root")
      #interaction[[count]] = append(interaction[[count]], as.character(tree_mat3$Current_node[tail(pos_new,1)]))
      #else
      #interaction[[count]] = append(interaction[[count]], root_of_age)
      pos_end=tail(pos_new,1)
      if(length(interaction[[count]])>1)
      {
        var_sampling_matrix_true<- rbind(var_sampling_matrix_true, tree_mat3[pos_end,])
        
        if(ind==1)
          interaction_summary<-rbind(interaction_summary,tree_mat3[pos_end,-c(1,ncol(tree_mat3))])
        else
          interaction_summary[ind,]<- rowSumBinary(interaction_summary[ind-1,],tree_mat3[pos_end,-c(1,ncol(tree_mat3))])
        var_sampling_matrix<- rbind(var_sampling_matrix, tree_mat3[pos_end,])
        var_sampling_matrix[vsm_row+ind,match(colnames(interaction_summary), colnames(var_sampling_matrix))]= interaction_summary[ind,]
        ind=ind+1
      }
    }
    vsm_row= nrow(var_sampling_matrix)
    count=count+1
  }
  #var_sampling_matrix <- var_sampling_matrix %>% distinct()  
  var_sampling_matrix <- var_sampling_matrix[-which(var_sampling_matrix[,1]=="root"),]
  return(list(interaction=interaction,var_sampling_matrix=var_sampling_matrix))
}

getInteractionMatrix <-function(output)
{
  ntree=length(output$nodes)
  var_num=ncol(output$data)-2
  var_names <- colnames(output$data)[-c(1:2)]
  Interaction_Score1 <- Interaction_Score2 <- var_interaction_score <- var_interaction_score_buf<- matrix(0, nrow=var_num, ncol=var_num)
  colnames(Interaction_Score1) <- rownames(Interaction_Score1) <- var_names
  colnames(Interaction_Score2) <- rownames(Interaction_Score2) <- var_names
  
  
  Age_Interaction_Score <- rep(0,var_num)
  names(Age_Interaction_Score) <- setdiff(var_names, "age")
  
  Interaction_Lists <- list()
  
  for(ind in 1:ntree)
  {
    
    #Interaction_Lists[[match(file, file_list)]][[ind]]<- list()
    #t<- party:::prettytree(output@ensemble[[ind]],names(output@data@get("input")))
    t <- output$nodes[[ind]]
    #plot(partykit::gettree(output, tree = ind))
    
    if(length(t$split)!=0)
    {
      if(names(t$info$p.value)!="age")
      {
        t_features <- GetPairVars(t,var_num, var_names)
        Var_Proportion <- Var_Interaction_Proportions(t_features[[1]])
        interaction <-  Var_Proportion[[1]]
        length(interaction)
        Interaction_Lists[[ind]] <- interaction
        
        var_sampling_matrix <- Var_Proportion[[2]]
        if(length(interaction)!=0)
        {
          var_interaction_mat <- matrix(0, nrow=length(interaction), ncol=var_num)
          var_sampling_mat<-  matrix(0, nrow=var_num, ncol=var_num)
          
          colnames(var_interaction_mat)=
            colnames(var_sampling_mat)=rownames(var_sampling_mat)=setdiff(var_names, "age")
          score_age <- rep(0,var_num)
          for( len in 1:length(interaction))
          {
            score <- rep(0,var_num)
            ind_match <- match(setdiff(interaction[[len]], "age"), names(Age_Interaction_Score))
            score[ind_match ]=1
            score_age <- score+score_age
            if(length(interaction[[len]])>2)
              var_interaction_mat[len,ind_match]=1
          }
          score_age[score_age>1]=1
          Age_Interaction_Score <- Age_Interaction_Score+score_age
          rownames(var_interaction_mat) <- paste("Interaction_Score_", 1:length(interaction), sep="")
          var_interaction_score <- var_interaction_score_buf <- t(var_interaction_mat)%*%var_interaction_mat
          
          #var_sampling_matrix<- as_tibble(var_sampling_matrix)
          var_sampling_matrix$Current_node <- as.character(var_sampling_matrix$Current_node)
          var_sampling_matrix[,setdiff(var_names,"age")] <- sapply(var_sampling_matrix[, setdiff(var_names,"age")], as.numeric)
          library(plyr)
          var_sampling_matrix <- ddply(var_sampling_matrix[,-match("Current_node", colnames(var_sampling_matrix))],"Immediate_parent",numcolwise(sum))
          row_substitute <- match(var_sampling_matrix$Immediate_parent, rownames(var_sampling_mat))
          var_sampling_mat[row_substitute,1:var_num]= as.matrix(var_sampling_matrix[,colnames(var_interaction_mat)])
          var_sampling_mat <- var_sampling_mat+t(var_sampling_mat)
          var_interaction_score <- var_interaction_score/var_sampling_mat
          var_interaction_score[is.na(var_interaction_score)]=0
          var_interaction_score[var_interaction_score==Inf]=0
        }
      }
    }
    Interaction_Score1 <- Interaction_Score1+var_interaction_score
    Interaction_Score2 <- Interaction_Score2+var_interaction_score_buf 
    
    
  }
  All_Interactions <- sapply(1:length(Interaction_Lists), function(tree)if(!is.null(Interaction_Lists[[tree]]))do.call(qpcR:::rbind.na, Interaction_Lists[[tree]]))
  All_Interactions <- do.call(qpcR:::rbind.na, All_Interactions)
  All_Interactions <- as.data.frame(All_Interactions)
  if(ncol(All_Interactions)==0)
    return(list( Interaction_Score1, All_Interactions))
  else
  {colnames(All_Interactions) <- paste("Node", 1:ncol(All_Interactions), sep="")
  All_Interactions <- All_Interactions%>% distinct()
  row_all_na <- sapply(1:nrow(All_Interactions), function(row)ifelse( all(is.na(All_Interactions[row,])), row,NA))
  if(length(na.omit(row_all_na))>0)
    All_Interactions <-All_Interactions[-na.omit(row_all_na),]
  
  All_Interactions_Stats <- cbind(All_Interactions, matrix(0, ncol=ntree, nrow=nrow(All_Interactions), dimnames = list(paste("Interaction_",1:nrow(All_Interactions), sep=""),c( paste("Tree_", 1:ntree, sep="")))))
  
  for(i in 1:length(Interaction_Lists))
  {
    occ <- sapply(1:nrow(All_Interactions), function(interaction) if(length(Interaction_Lists[[i]])!=0){sum(sapply(1:length(Interaction_Lists[[i]]), function(tree)all(na.omit(unlist(All_Interactions[interaction,])) %in% Interaction_Lists[[i]][[tree]])))})
    if(!is.null(unlist(occ)))
      All_Interactions_Stats[,paste("Tree_", i, sep="")] =occ
  }
  
  
  All_Interactions_Stats[ is.na( All_Interactions_Stats)]=""
  All_Interactions_Stats[,paste("Tree_", 1:ntree, sep="")][ which( All_Interactions_Stats[, paste("Tree_", 1:ntree, sep="")]>1, arr.ind=TRUE)]=1
  All_Interactions_Stats$Forest_Score <- rowSums( All_Interactions_Stats[, paste( "Tree_",1:ntree, sep="")])
  All_Interactions_Stats <- All_Interactions_Stats[ order( All_Interactions_Stats$Forest_Score),]
  All_Interactions_Stats <- All_Interactions_Stats[,-match(paste("Tree_",1:ntree,sep=""), colnames(All_Interactions_Stats))]
  
  
  return( list( Interaction_Score1, All_Interactions_Stats))
  }
}


RowContain <- function(set, mat) sapply(1:nrow(mat), function(r)ifelse(all(set %in% mat[r,]), r, NA))

GenerateInteractionList <- function(Interaction_List, Importance_Score_list, counter)
{
  
  All_Interactions_Stats <- do.call(qpcR:::rbind.na, Interaction_List)
  All_Interactions_Stats <- All_Interactions_Stats [,-1]
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
    mat_to_check <- mat_to_check[,-1]
    mat_to_check <- mat_to_check %>% filter(Forest_Score>0)
    mat_to_check[mat_to_check==""]="z"
    mat_to_check[,-ncol(mat_to_check)] <- t(apply( mat_to_check[,-ncol(mat_to_check)], 1, grr::sort2))
    mat_to_check[mat_to_check=="z"]=""
    
    if(length(which(colSums(mat_to_check=="")==nrow(mat_to_check)))>0)
      mat_to_check <- mat_to_check[, which(colSums(mat_to_check=="")==nrow(mat_to_check))]
    
    system.time(pos_match <- matchRow(nrow(mat_to_check), as.matrix(All_Interactions_Stats[,1:(ncol(mat_to_check)-1)]), as.matrix(mat_to_check[, -ncol(mat_to_check)])))
    #system.time(pos_match <-vapply(1:nrow(mat_to_check), function(row)which(compareToRow(as.matrix(All_Interactions_Stats[,1:(ncol(mat_to_check)-1)]),as.vector(unlist(mat_to_check[row,-ncol(mat_to_check)])))), numeric(1)))
    
    All_Interactions_Stats[pos_match,paste("Forest_",i,sep="")] = mat_to_check$Forest_Score
    #All_Interactions_Stats[[paste("Forest_",i,sep="")]][All_Interactions_Stats[[paste("Forest_",i,sep="")]]>1]=1
  })
  
  All_Interactions_Stats$Sum_Forest_Score <- sapply( 1:nrow(All_Interactions_Stats), function(i) sum( as.numeric(All_Interactions_Stats[ i, paste("Forest_", 1:NiterChild,sep="")]), na.rm=TRUE) )
  All_Interactions_Stats$Median_Forest_Score <- rowMedians(as.matrix(All_Interactions_Stats[,paste("Forest",1:NiterChild, sep="_")]))
  All_Interactions_Stats <- All_Interactions_Stats%>% distinct()
  All_Interactions_Stats <- All_Interactions_Stats[order(All_Interactions_Stats$Median_Forest_Score),]
  All_Interactions_Stats <- All_Interactions_Stats
  #All_Interactions_Stats<-All_Interactions_Stats[which(All_Interactions_Stats$Median_Forest_Score>2),]
  return(All_Interactions_Stats)
}

ComputeMedian <- function(temp)
{a <- do.call(abind, c(temp, list(along=3)))
Interaction_Score <- apply(a, 1:2, function(m)median(m, na.rm=TRUE))

if(max(Interaction_Score)==0)
  return(NA)
else
{i=1 
while(i < nrow(Interaction_Score))
{
  if(all(Interaction_Score[i,]==0))
  {Interaction_Score<-Interaction_Score[-i,-i]
  i=i-1
  }
  i=i+1
}
Interaction_Score <- Interaction_Score[-ncol(Interaction_Score), -ncol(Interaction_Score)]
return(Interaction_Score)
}
}

null <- function(n)ifelse(length(n)==0,NA,n)
addq <- function(x) paste0("`", x, "`")

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

Max.Test <- function(Cluster)
{
  
   #For each combination we need to find the beta that is the effect due to the age or the temporal effect
  # Y~\betaxt + \alpha + \epsilon
  SNP_Level_Summary <- table(Cluster$age, Cluster$SNP_Clubbed)
  To_Remove <- c(names(which(SNP_Level_Summary[1,]==0)), names(which(SNP_Level_Summary[2,]==0)))
  if(length(To_Remove)>=1)
    Cluster <- Cluster %>% filter(!(SNP_Clubbed%in%To_Remove))
  Cluster$SNP_Clubbed <- factor(Cluster$SNP_Clubbed)
  SNP_levels <- levels(Cluster$SNP_Clubbed)
  reference <- SNP_levels[1]
  Beta_Estimate <-Beta_SE <- rep(0, length(SNP_levels))
  T <- Stderror<-  rep(0, (nlevels(Cluster$SNP_Clubbed)-1))
  df <- rep(0, length(SNP_levels))
  Age_levels <- levels(Cluster$age)
  
  for(snp_level in SNP_levels)
  {
    i <- match(snp_level, SNP_levels)
    #Filtering the data have a particular genotype combination
    Cluster_reference <- Cluster %>% filter(SNP_Clubbed==reference)
    Cluster_snp_level <- Cluster %>% filter(SNP_Clubbed==snp_level)
    #Getting the number of observations across each age in the reference genotype combination
    N0 <- table(Cluster %>% filter(SNP_Clubbed==reference)%>%dplyr::select("age"))
    #Getting the number of observations across each age in the ith genotype combination
    N <- table(Cluster_snp_level$age)
    #Fitting the linear regression model
    fit_level <- lm(value~age, data=Cluster_snp_level)
    #Estimating the longitudinal effect on the genotype combination i
    Beta_Estimate[i]=coef(fit_level)[2]
    #
    #sigma12 <- sum(var(Cluster%>% filter(SNP_Clubbed==reference&age==Age_levels[1])%>%select("value"))/N0[1],
    #    var(Cluster%>% filter(SNP_Clubbed==reference&age==Age_levels[2])%>%select("value"))/N0[2])
    #Estimating the variance of the longitudinal effect on the genotype combination i
    s1 <- sqrt(var(  Cluster_snp_level$value[which(Cluster_snp_level$age==Age_levels[1])]))
    s2 <- sqrt(var(  Cluster_snp_level$value[which(Cluster_snp_level$age==Age_levels[2])]))
    Beta_SE[i] <-sqrt(s1^2/N[1]+s2^2/N[2])
    
    #t_reference <- t.test(value~age, Cluster_reference)
    #s1 <- sqrt(var(  Cluster_reference$value[which(Cluster_reference$age==1)]))
    #s2 <- sqrt(var(  Cluster_reference$value[which(Cluster_reference$age==2)]))
    #Beta_SE[ match(reference, SNP_levels)] <-sqrt(s1^2/N0[1]+s2^2/N0[2])
    #Beta_SE[ match(reference, SNP_levels)] <- t_reference$stderr
    
    
    DF_num <- (s1^2/N[1]+s2^2/N[2])^2
    DF_den <- (s1^2/N[1])^2/(N[1]-1)+(s2^2/N[2])^2/(N[2]-1)
    df[i] <- DF_num/DF_den
    if(snp_level!=SNP_levels[1])
    {
      Stderror[i-1] <- sqrt(Beta_SE[i]^2+Beta_SE[match(reference, SNP_levels)]^2)
      T[ i-1]=(Beta_Estimate[i]-Beta_Estimate[match(reference, SNP_levels)])/Stderror[i-1]
    }
    
  }
  
  Ti.Table <- data.frame(diff_combi=Beta_Estimate[-1]-Beta_Estimate[1], Stderror= Stderror, estimate=T, DF=df[-1]+df[1] )
  rownames( Ti.Table) <- setdiff(SNP_levels, reference)
  Ti.Table <- Ti.Table%>%mutate(pval.unadj= 1-(pt(abs(T), DF)-pt(-abs(T), DF)))
  Ti.Table <- Ti.Table%>%mutate(pval.adj= 1-(pt(max(abs(T)), DF)-pt(-max(abs(T)), DF)))
  R=diag(nlevels(Cluster$SNP_Clubbed)-1) 
  R[upper.tri(R) | lower.tri(R)] = 0.5
  S1 <- var(Cluster%>%filter(age==Age_levels[1])%>%dplyr::select("value"))
  S2 <- var(Cluster%>%filter(age==Age_levels[2])%>%dplyr::select("value"))
  N1 <- table(Cluster$age)[1]
  N2 <- table(Cluster$age)[2]
  df_num <- (S1/N1+S2/N2)^2
  df_den <- ((S1/N1)^2/(N1-1))+((S2/N2)^2/(N2-1))
  df_final <- sqrt(2)*df_num/df_den  
  #df_final <- sum(Ti.Table$DF)/(nlevels(Cluster$SNP_Clubbed)-1) 
  pval <- 1-pmvt(lower=rep(-max( abs(Ti.Table$estimate)),length(SNP_levels)-1), upper=rep(max( abs(Ti.Table$estimate)),length(SNP_levels)-1), delta=0, df=round(df_final,0), corr=R)
  pval2 <-1- pmvnorm(lower=rep(-max( abs(Ti.Table$estimate)),length(SNP_levels)-1), upper=rep(max(abs(Ti.Table$estimate)),length(SNP_levels)-1), mean=rep(0, length(SNP_levels)-1)
          , sigma=R, algorithm = GenzBretz())
  return(list(Ti.Table, pval, pval2))
}


Final.Contrast.MaxTest2 <- function(listCluster,Cluster1, niter){
  #random_cluster_pval <- sapply(1:niter, function(i)median(sapply( listCluster[[i]], function(clustsampledata)compute_pval(clustsampledata,filtercomb,contrast,slope,eqVar,g, Data_Phenotype))))
  system.time({
    n_cores <- detectCores(logical=FALSE)
    cl <- makeCluster(n_cores-1)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library("dplyr"))
    clusterEvalQ(cl, library("plyr"))
    clusterEvalQ(cl, library("mvtnorm"))
    clusterExport(cl, c("Max.Test2","listCluster", "niter", "pmvt"))
    random_cluster_pval <- parLapply(cl , 1:niter,  function(i)median(sapply( listCluster[[i]], function(clustsampledata){p <- Max.Test2(clustsampledata);return(p)})))
  })
  #stopCluster(cl)
  # If NaNs are produced 
  random_cluster_pval <- unlist(random_cluster_pval[!is.na(random_cluster_pval)])
  
 desired_cluster_pval <- median(sapply(1:length(Cluster1), function(ind) 
  {p <-Max.Test2(Cluster1[[ind]]);return(p)}))
  
  print(desired_cluster_pval)
  prop <- sum(random_cluster_pval<desired_cluster_pval)/length(random_cluster_pval)
  return(prop)
}

Final.Contrast.MaxTest <- function(listCluster,Cluster1, niter){
  #random_cluster_pval <- sapply(1:niter, function(i)median(sapply( listCluster[[i]], function(clustsampledata)compute_pval(clustsampledata,filtercomb,contrast,slope,eqVar,g, Data_Phenotype))))
  system.time({
    n_cores <- detectCores(logical=FALSE)
    cl <- makeCluster(n_cores-1)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library("dplyr"))
    clusterEvalQ(cl, library("plyr"))
    clusterEvalQ(cl, library("mvtnorm"))
    clusterExport(cl, c("Max.Test","listCluster", "niter", "pmvt"))
    random_cluster_pval <- parLapply(cl , 1:niter,  function(i)sapply( listCluster[[i]], function(clustsampledata){p <- Max.Test(clustsampledata);return(p)}))
    random_cluster_pval1 <- sapply(1:niter, function(i)median(unlist(random_cluster_pval[[i]][2,])))
    random_cluster_pval2 <- sapply(1:niter, function(i)median(unlist(random_cluster_pval[[i]][3,])))
    
  })
  
 
  #stopCluster(cl)
  # If NaNs are produced 
  #random_cluster_pval <- unlist(random_cluster_pval[!is.na(random_cluster_pval)])
  
  desired_cluster_pval <- sapply(1:length(Cluster1), function(ind) 
  {p <-Max.Test(Cluster1[[ind]]);return(p)})
  desired_cluster_pval1 <- median(unlist(desired_cluster_pval[2,])) 
  desired_cluster_pval2 <- median(unlist(desired_cluster_pval[3,])) 
  
  #print(desired_cluster_pval)
  prop1 <- sum(random_cluster_pval1<desired_cluster_pval1)/length(random_cluster_pval1)
  prop2 <- sum(random_cluster_pval2<desired_cluster_pval2)/length(random_cluster_pval2)
  return(c(prop1, prop2))
}

plotSNPAgingInteraction <- function(data_epistasis, snp_list)
{
  
  data_epistasis$Genotype_Combination <- factor(sapply(1:nrow(data_epistasis),FUN=function(x)paste(unfactor(data_epistasis[x,snp_list]),collapse = "")))
  g <- list()
  xlabs <- paste(levels(data_epistasis$Genotype_Combination),"\n (N1=",table(data_epistasis$age, data_epistasis$Genotype_Combination)[1,],")","\n (N2=",table(data_epistasis$age, data_epistasis$Genotype_Combination)[2,],")",sep="")
  g[[1]] <- ggplot(data_epistasis, aes(x=Genotype_Combination, y=value, fill=age))+scale_fill_viridis(discrete = TRUE, alpha=0.6)+scale_x_discrete(labels=xlabs)+
     geom_boxplot()+theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size=11))
  for(i in 1:length(snp_list))
  {
    j=i+1
    xlabs <- paste(levels(data_epistasis[[snp_list[i]]]),"\n (N1=",table(data_epistasis$age, data_epistasis[[snp_list[i]]])[1,],")","\n (N2=",table(data_epistasis$age, data_epistasis[[snp_list[i]]])[2,],")",sep="")
    g[[i+1]] <- ggplot(data_epistasis, aes(x=.data[[snp_list[i]]], y=value, fill=age))+scale_fill_viridis( discrete = TRUE, alpha=0.1)+scale_x_discrete(labels=xlabs) +
      geom_boxplot()+theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(size=11))+xlab(snp_list[i])
  }
  fit1 <- lm(as.formula(paste("value~age*",paste(addq(snp_list), collapse="+"))), data= data_epistasis)
  fit2 <- lm(as.formula(paste("value~","age*Genotype_Combination", sep="")), data= data_epistasis)
  an <- anova(fit1, fit2, test="Chisq")
  lay <- rbind(2:(length(snp_list)+1),rep(1, length(snp_list)))
  #g[[i+2]]<- grid.table(an)%>%kable()
  grid.arrange(grobs=g, layout_matrix = lay,  top = textGrob(paste("Testing the effect of the n-way interaction between", paste(snp_list, collapse = ","), sep=""),gp=gpar(fontsize=10,font=2)))
  
  print(an)
}


Max.Test2 <- function(Cluster)
{
  SNP_Level_Summary <- table(Cluster$age, Cluster$SNP_Clubbed)
  To_Remove <- c(names(which(SNP_Level_Summary[1,]==0)), names(which(SNP_Level_Summary[2,]==0)))
  if(length(To_Remove)>=1)
    Cluster <- Cluster %>% filter(!(SNP_Clubbed%in%To_Remove))
  Cluster$SNP_Clubbed <- factor(Cluster$SNP_Clubbed)
  SNP_levels <- levels(Cluster$SNP_Clubbed)
  reference <- levels(SNP_levels)[1]
  #For each combination we need to find the beta that is the effect due to the age or the temporal effect
  # Y~\betaxt + \alpha + \epsilon
  T <-  estimate <- rep(0, length(SNP_levels))
  df <- pval <- rep(0, length(SNP_levels))
  #t_full <- t.test(value~age, Cluster)
  #mu <- diff(t_full$estimate)
  Cluster$value <- Cluster$value-predict(lm(value~age, data=Cluster))
  for(snp_level in SNP_levels)
  {
    i <- match(snp_level, SNP_levels)
    #Filtering the data have a particular genotype combination
    #Cluster_reference <- Cluster %>% filter(SNP_Clubbed==reference)
    Cluster_snp_level <- Cluster %>% filter(SNP_Clubbed==snp_level)
    Cluster_snp_level_t <- t.test(value~age, data=   Cluster_snp_level,mu=0, alternative="two.sided", var.equal=FALSE)
    T[i]=Cluster_snp_level_t$statistic
    pval[i]=Cluster_snp_level_t$p.value
    df[i]=Cluster_snp_level_t$parameter
    estimate[i]=diff(Cluster_snp_level_t$estimate)
  }
    tmax <- max(abs(T))
    pval <- 1- prod(sapply(df, function(d)pt(tmax,df=d)-pt(-tmax,df=d)))
    return(pval)
}


Sample_within_cluster_ForRandomSamples <- function(snp_sample, Phenotype_Data, SNP_Sampling_List)
{
  snp_sample <- sample(setdiff(SNP_Sampling_List, snp_sample),size = length(snp_sample), replace = F)
  
  lst_Random_SNP_clusters <- Phenotype_Data[, c("value", "age", "strain", unlist(  snp_sample))]
  lst_Random_SNP_clusters$SNP_Clubbed <- sapply(1:nrow(lst_Random_SNP_clusters),function(x)paste(lst_Random_SNP_clusters[x,4:ncol(lst_Random_SNP_clusters)],collapse = ""))
  lst_Random_SNP_clusters$age <- factor( lst_Random_SNP_clusters$age)
  Seeds=sample(1:100,10)
  data <- lst_Random_SNP_clusters
  lst_Random_SNP_clusters_Sampled <-list()
  #test1 <- list()
  system.time(for(seed in Seeds )
  {
    s = match(seed, Seeds)
    lst_Random_SNP_clusters_Sampled[[s]] <- Sample_within_cluster(lst_Random_SNP_clusters, filtercomb = "no", g=30, seed)
    lst_Random_SNP_clusters_Sampled[[s]]$SNP_Clubbed <- factor(lst_Random_SNP_clusters_Sampled[[s]]$SNP_Clubbed)
    lst_Random_SNP_clusters_Sampled[[s]]$age <- factor(lst_Random_SNP_clusters_Sampled[[s]]$age)
    #test1[[match(seed, Seeds)]] <- compute_pval2( Data_Phenotype_Cluster_Sampled,filtercomb = "no",contrast="second",slope="yes", eqVar="no", g=1, Phenotype_Data)
    #test1[[match(seed, Seeds)]][,1:5] <- round( test1[[match(seed, Seeds)]][ ,1:5], 2)
    #test1[[match(seed, Seeds)]]$p.value.unadj <- format.pval(test1[[match(seed, Seeds)]]$p.value.unadj)
    #test1[[match(seed, Seeds)]]$p.value <- format.pval(test1[[match(seed, Seeds)]]$p.value)
    #print(kable(test1[[match(seed, Seeds)]]) %>%  kable_styling(bootstrap_options = c("striped", full_width = F)))
  })
  return(lst_Random_SNP_clusters_Sampled)
}

Sample_within_cluster_ForRandomSamples_woBootstrap <- function(snp_sample, Phenotype_Data, SNP_Sampling_List)
{
  snp_sample <- sample(setdiff(SNP_Sampling_List, snp_sample),size = length(snp_sample), replace = F)
  
  lst_Random_SNP_clusters <- Phenotype_Data[, c("value", "age", "strain", unlist(  snp_sample))]
  lst_Random_SNP_clusters$SNP_Clubbed <- sapply(1:nrow(lst_Random_SNP_clusters),function(x)paste(lst_Random_SNP_clusters[x,4:ncol(lst_Random_SNP_clusters)],collapse = ""))
  lst_Random_SNP_clusters$age <- factor( lst_Random_SNP_clusters$age)
  return(lst_Random_SNP_clusters)
}


## Function that filters out combinations that are not present in both ages ##g 
## I added g as a parameter as here one can easily mention the number of observations threshold
## one want while building the test statistics. This can help you avoid calling a separate function
## RemoveObs
FilterCluster <- function( Cluster, g){
  
  Cluster$age <- factor( Cluster$age)
  Cluster$SNP_Clubbed <- factor( Cluster$SNP_Clubbed)
  # Frequency of SNP_combinations
  table_stat <- table( Cluster$SNP_Clubbed, Cluster$age)
  # Filter out combinations that are not present in both ages
  if( length( names( which(table_stat[,1]<=g | table_stat[,2]<=g)))!= 0){
    Cluster <- Cluster %>% filter( !(SNP_Clubbed %in% names( which( table_stat[,1]<=g | table_stat[,2]<=g))))
    Cluster$SNP_Clubbed <- factor( Cluster$SNP_Clubbed)
  }
  return( Cluster)
  
}

Sample_within_cluster <- function( Cluster, filtercomb, g,  Seed){
  
  x <- Cluster # It will be needed in an if statement below
  #Data_Phenotype <- Data_Phenotype ##reduced data
  T.Statistics <- list()
  
  # If we want to filter out 1's combns 
  if (filtercomb == "yes"){
    # Filter Cluster
    Cluster <- FilterCluster(Cluster, g)
    # Filter 1's combination
    Cluster <- FilterComb(Cluster)
  }
  else
    Cluster <- FilterCluster(Cluster, g)
  
  # Frequency of SNP combinations in age 1 and age 2
  table_stat_orig <- table(Cluster$SNP_Clubbed , Cluster$age)
  combns <- rownames( table_stat_orig)
  Cluster_Sampled <- list()
  for(comb in combns)
  {
    c <- match(comb, combns)
    Cluster_New <- Cluster %>% filter(SNP_Clubbed==comb)
    if(table_stat_orig[c,1]<80 | table_stat_orig[c,2] <80)
      Cluster_Sampled[[c]]=Cluster_New
    else
    {
      Cluster_Comb <- Cluster_New  
      Cluster_Comb_Age1 <- Cluster_Comb %>% filter(age==levels(Cluster_Comb$age)[1])
      Cluster_Comb_Age2 <- Cluster_Comb %>% filter(age==levels(Cluster_Comb$age)[2])
      Cluster_Sampled[[c]]=rbind(Cluster_Comb_Age1[sample(nrow(Cluster_Comb_Age1),size=80,replace=FALSE),], Cluster_Comb_Age2[sample(nrow(Cluster_Comb_Age2),size=80,replace=FALSE),])
    }
  }
  Cluster_Sampled <- do.call("rbind", Cluster_Sampled)
  return(Cluster_Sampled)
}