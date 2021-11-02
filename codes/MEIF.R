
############################################################################################################
################################The creators of MEIF are as follows########################################
#############################													                        ##############################
###########################		    		Saswati Saha		 	                       #############################
##########################	   TAGC, Inserm, AMU, Marseille			  	              ##########################


###	Revised: 02 Nov 2021	###
##############################################################################################################################

##########################I have just modified the merf code provided by the original users (Hajjem et al.,),#################
#################           such that it can incorporate a conditional forest instead of                      ################ 
###############              random forest. This is needed to avoid missing value                              ###############
##############                 treatment for SNP data. Also the following code                                  ##############
#############                     can contain 3 components:a cforest                                              ############
###########                       with fixed effects, random effects,                                              ###########
###########                       fixed effects outside the cforest,                                                ##########
##########                         parallelization of the  cforest.                                                   ########
##########                                 .                                                                          ########

#############################################################################################################
###	Description of MEIF arguments	###
###########################################

#xnam		#A charcter vector of p columns, corresponding to the names of the p fixed effect covariates (as they appear in the learning dataset).

#MERF.lDB 	#The learning dataset: a dataframe of N rows (i.e. N level-one observations) and (p+2) columns, 
#where the column "cluster.id" corresponds to the variable uniquely identifying the n clusters (or level-two obervations), 
#the column "Y" corresponds to the name of the continuous response variable,
#and the other p columns correspond to the p fixed effects covariates.	

#ni			#A vector of n columns, where ni[i] corresponds to the size of cluster i, 
#for i = 1, ..., n (ATTENTION: should keep the same order as in MERF.lDB).

#Zi			#A list of n matrices Zi[[i]] of dimension (ni[i] X q), where  
#Zi[[i]] corresponds to the q random effects covariates values of the ni[i] observations nested within cluster i, for i= 1, ..., n. 
#Note that q=1 when there is only a random intercept, i.e. no random effect covariates.

#Xi     #A list of n matrices Zi[[i]] of dimention(ni[[i]] x nx) where
#Xi[[i]] corresponds to the fixed effects covariates that we donot wish to include in the random forest for cluster i where i =1, ..., n.

#Yi			#A list of n vectors Yi[i] of ni[i] rows, where 
#Yi[i] corresponds to the response values of the ni[i] observations nested within cluster i, for i= 1, ..., n.

#ntree			#ntree argument of randomForest function, i.e., number of trees to grow (e.g. 300).

#mtry			#mtry argument of randomForest function, i.e., number of variables randomly sampled as candidates at each split (e.g. 3).

#nodesize		#nodesize argument of randomForest function, i.e., minimum size of terminal nodes (e.g.5).

#sigmasqzero = NULL	#Starting value of s2, where the covariance matrix of the errors Ri = s2Ini , for i = 1, ..., n. 
#Default: if( is.null(sigmasqzero) ) sigmasqzero <- 1.

#Dzero = NULL		#Starting value of the covariance matrix of the random effects. 
#Default:
#if( is.null(Dzero) ){
#	Dzero <- diag(0.01,nrow=q,ncol=q)
#}

#bizero = NULL		#Starting values of the random effects: a list of n matrices of q × 1 unknown vector of random effects. 
#Default:
#if( is.null(bizero) ){
#bizero <- list(); length(bizero)<- n
#	for(i in 1:n) bizero[[i]] <- matrix(0,nrow=q,ncol=1)}
#}

#cizero = NULL		#Starting values of the fixed effects( not part of random forest): a list of n matrices of nx × 1 unknown vector of fixed effects. 
#Default:
#if( is.null(cizero) ){
#cizero <- list(); length(cizero)<- n
#	for(i in 1:n) cizero[[i]] <- matrix(0,nrow=nx,ncol=1)}
#}
#F.niter	#The number of iterations forced to avoid early stopping (e.g. 100).

#max.niter	#Maximum number of iterations, in addition to the F.niter forced iterations (e.g. 300).

#smallest.Jump.allowed 	#A given small value (e.g. 1e-4). 

#verbose	#Logical. Should R report extra information on progress?
#threads	#number. Number of threads used for parallelization



####################################	BEGIN THE SCRIPT		#####################################

#setwd("/home/...")
#sink("Routput.doc")
set.seed(321)

########################################################################################################################
###############################                   ### MERF function ###  	         ##############################
##################		fits a mixed effects random forest of regression trees model		###############
########################################################################################################################


##############################
# 0 - Load librairies
##############################
library( randomForest)
library( lmerTest)
library( parallel)
library( dplyr)
library( partykit)
library (MASS)
library(snowfall)
library(partykit)
############################## 
# 1 - Source file 
##############################


cforestmt<-function(formula, data = list(), subset = NULL, weights = NULL,weight_variable=NULL, control = ctree_control(),mtry=mtry, ntree=ntree, xtrafo = NULL, ytrafo = NULL, scores = NULL, threads=6) {
  
  if(ntree<threads) {    # if there are less trees than threads single thread
    return(cforest_new(formula, data = data, subset=subset,ntree=ntree, mtry=mtry, weights=weights,weight_variable=weight_variable, control=control, xtrafo=xtrafo, ytrafo=ytrafo, scores=scores))
  }
  
  # round off threads
  fsize=ntree/threads
  if(fsize-round(fsize)!=0) {
    fsize=ceiling(fsize)
    message("Rounding forest size to ",fsize*threads)
  }
  ntree=as.integer(fsize)
  
  # run forests in parallel
  sfInit(parallel=T, cpus=threads, type="SOCK")
  sfClusterEval(library(partykit))
  sfExport('formula','data','subset','weights','weight_variable','control','xtrafo','ytrafo','scores', 'mtry', 'ntree', "cforest_new", ".start_subset", "constparties", "ctree_new", ".y2infl")
  fr<-sfClusterEval(cforest_new(formula, data = data, subset=subset, weights=weights,weight_variable=weight_variable, mtry=mtry,ntree=ntree, control=control,  ytrafo=ytrafo, scores=scores))
  sfStop()
  
  # combine/append forest
  fr[[1]]$nodes<-unlist(lapply(fr,function(y) {y$nodes}),recursive=F)
  fr[[1]]$info<-unlist(lapply(fr,function(y) {y$info}),recursive=F)
  fr[[1]]$weights<-unlist(lapply(fr,function(y) {y$weights}),recursive=F)
  
  #first forest has result
  return(fr[[1]])
}


MERF  <- function( xnam,
                   MERF.lDB, 
                   ni,
                   Zi,
                   Yi, 
                   Xi,
                   ntree,
                   mtry,
                   nodesize,
                   sigmasqzero = NULL,
                   Dzero = NULL,
                   bizero = NULL,
                   cizero = NULL,
                   weight_variable=NULL,
                   F.niter,
                   max.niter,
                   smallest.Jump.allowed,
                   threads=6,
                   verbose = TRUE)
{
  #Parameters values 
  n <- length(ni)
  N <- sum(ni)
  q <- dim(Zi[[1]])[2] 	# q=1 in random intercept case
  nx <- dim(Xi[[1]])[2]  
  
  #Initial values of sigmasqzero, Dzero, and bizero
  if( is.null(sigmasqzero) ) sigmasqzero <- 1
  else sigmasqzero <- sigmasqzero
  
  if( is.null(Dzero) ){
    Dzero <- diag( 0.01 , nrow = q, ncol = q)
  }
  else Dzero <- Dzero
  
  if( is.null(bizero) ){
    bizero <- list(); length(bizero) <- n
    for(i in 1:n) bizero[[i]] <- matrix( 0, nrow = q, ncol = 1)
  }	
  else bizero <- bizero
  
  if( is.null(cizero) ){
    cizero <- matrix(0, nrow = nx, ncol = 1)}
  else cizero <- cizero
  
  #iter number
  r <- 1
  if (verbose) 	
    message("MERF iter no: ", r)
  
  #transformed outcome, star.Yi[[r]][[i]], initialized with the original values
  star.Yi <- list() 
  
  for(i in 1:n)
    star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bizero[[i]] - Xi[[i]] %*% cizero
  
  MERF.lDB$star.Yi <- rep(0, nrow(MERF.lDB))
  
  for(i in 1:n)
    MERF.lDB[which(MERF.lDB$cluster.id==levels(MERF.lDB$cluster.id)[i]),]$star.Yi <- star.Yi[[i]]
  
  
  
  #one STD random forest
  ######################
  
  #MERF.lDB$star.Yi <- unlist(star.Yi) 
  rm(star.Yi) ;  gc(verbose=FALSE)
  addq <- function(x) paste0("`", x, "`")
  fit.rf.formula <-  as.formula( paste(" star.Yi ~ ", paste( addq(colnames( xnam)), collapse = "+")))
  
  #MERF.lDB$f.pred <- rep(0, nrow(MERF.lDB))
  ##fitting a cforest instead of random forest
  fit.rf <- cforestmt( formula = fit.rf.formula , data = MERF.lDB	, weight_variable=weight_variable, mtry = mtry, ntree = ntree, control = ctree_control(teststat="quad", testtype="Univ", mincriterion=0.95, minsplit = 50, minbucket = 50, maxsurrogate = min(1, ncol( xnam))), threads = threads)
  
  
  MERF.lDB$f.pred  <- fit.rf$fitted$`(response)`
  #in matrix format
  fixed.pred <- list()	
  fixed.pred <- split(MERF.lDB, MERF.lDB$cluster.id) 	
  for(i in 1:n)fixed.pred[[i]] <- as.matrix( subset( fixed.pred[[i]] ,select = f.pred), ncol=1)
  
  #random part
  ############
  #random effects parameters in list format
  bi <- list(list()) ; length(bi) <- r
  for(i in 1:n) bi[[r]][[i]] <- bizero[[i]] 	
  #print("bizero"); print(bizero)
  rm(bizero) ; gc(verbose=FALSE)
  
  #fixed part
  ############
  #random effects parameters in list format
  ci <- list() ; length(ci) <- r
  ci[[r]] <- cizero	
  #print("bizero"); print(bizero)
  rm(cizero) ; gc(verbose=FALSE)
  
  #level-1 variance component
  #residuals
  epsili <- list()
  for(i in 1:n)
    epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]] - Xi[[i]] %*% ci[[r]]
  
  sigma.sq <- vector( mode = "numeric") ; length( sigma.sq ) <- r
  sigma.sq[r] <- sigmasqzero
  #print("sigmasqzero") ; print(sigmasqzero)			
  rm(sigmasqzero) ; gc(verbose=FALSE)
  #message("sigmasq of current micro iter", sigma.sq[r] )
  
  #level-2 variance component
  D <- list() ; length(D) <- r
  D[[r]] <- Dzero		#!!!Dzero <- diag(x=0.01, nrow=q, ncol = q)
  #print("Dzero") ; print(Dzero)	
  rm(Dzero) ; gc(verbose=FALSE)
  #message("D of current micro iter: ", D[[r]] )
  
  #level-1 and level-2 variance components (or typical or total variance)
  Vi <- list() 
  inv.Vi <- list(list()) ; length(inv.Vi) <- r
  
  for(i in 1:n){
    Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]]) + sigma.sq[r]*diag( x = 1, nrow = ni[i], ncol = ni[i])
    if(q==1)
      inv.Vi[[r]][[i]] <- 
        (1/sigma.sq[r]) * ( diag( rep( 1, ni[i]))
                            -(( as.numeric( D[[r]])/sigma.sq[r])/( 1+ni[i]*( as.numeric( D[[r]])/sigma.sq[r])))*matrix( rep( 1,(ni[i])^2), ncol=ni[i], nrow=ni[i]) )
    else inv.Vi[[r]][[i]] <- solve( Vi[[i]])
  }
  
  Vi <- list(NULL) 
  #inv.Vi[[r-1]] <- list(NULL) #not to run at step 0
  
  #the generalized log-likelihood (GLL) 
  GLL <- vector(mode="numeric") ; length(GLL) <- r
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i] <- t( epsili[[i]]) %*% solve( sigma.sq[r]*diag( x = 1, nrow = ni[i], ncol = ni[i])) %*% epsili[[i]]+
    t( bi[[r]][[i]]) %*% solve( D[[r]]) %*% bi[[r]][[i]]+
    log( det(D[[r]]))+
    log( det( sigma.sq[r]*diag( x = 1, nrow = ni[i], ncol = ni[i])))
  GLL[r] <- sum(term)
  rm(term)
  gc(verbose=FALSE)
  
  #convergence criterion
  Jump <- rep(NA, r) 		#at this first iteration Jump = NA
  convergence.iter <- rep(NA, r) 	#at this first convergence.iter = NA
  
  ####################
  ####	STEP 1	####        
  ####################
  
  #update iteration number r
  r <- r+1
  if (verbose) 	
    message("MERF iter no: ", r)
  
  #update the length of the different lists
  
  length(sigma.sq) <- r
  length(D) <- r
  length(inv.Vi) <- r
  length(bi) <- r
  length(ci) <- r
  length(GLL) <- r
  
  length(Jump) <- r
  length(convergence.iter) <- r
  
  #update the transformed outcome, star.Yi
  star.Yi <- list() 
  for(i in 1:n){
    star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	- Xi[[i]] %*% ci[[r-1]]
  }	
  
  #one STD random forest
  ######################
  
  for(i in 1:n)
    MERF.lDB[which(MERF.lDB$cluster.id==levels(MERF.lDB$cluster.id)[i]),]$star.Yi <- star.Yi[[i]]
  rm(star.Yi) ;  gc(verbose=FALSE)
  

  
  fit.rf <- cforestmt( formula = fit.rf.formula , data = MERF.lDB	,weight_variable=weight_variable,  mtry = mtry, ntree = ntree, control = ctree_control(teststat="quad", testtype="Univ", mincriterion=0.95, minsplit = 50, minbucket = 50, maxsurrogate = min(1, ncol( xnam))), threads = threads)
  
  
  MERF.lDB$f.pred  <-  fit.rf$fitted$`(response)`
  
  
  
  #in matrix format (from random forest)
  fixed.pred <- list()	
  fixed.pred <- split( MERF.lDB, MERF.lDB$cluster.id) 	
  for(i in 1:n)fixed.pred[[i]] <- as.matrix( subset(fixed.pred[[i]] , select=f.pred), ncol=1)
  
  #fixed part
  ###########
  A = t(Xi[[1]]) %*% inv.Vi[[r-1]][[1]] %*% Xi[[1]]
  B = t(Xi[[1]]) %*% inv.Vi[[r-1]][[1]] %*% ( Yi[[1]] - fixed.pred[[1]])
  for(i in 11:20)
    { A= A + t(Xi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% Xi[[i]]
      B=B + t(Xi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% ( Yi[[i]] - fixed.pred[[i]])
    }
  ci[[r]] <- ginv(A) %*% B
  #random	part
  ############
  for(i in 1:n)
    bi[[r]][[i]] <- D[[r-1]] %*% t( Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% ( Yi[[i]] - fixed.pred[[i]] - Xi[[i]] %*% ci[[r]])
  
  bi[r-1] <- list(NULL)		 
  
  ####################
  ####	STEP 2	####         
  ####################
  
  #level-1 variance component
  #residuals
  epsili <- list()
  for(i in 1:n)
    epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]] - Xi[[i]] %*% ci[[r]]
  
  
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i] <- crossprod(epsili[[i]]) + 
    sigma.sq[r-1] * ( ni[i] - sigma.sq[r-1]* sum( diag( inv.Vi[[r-1]][[i]])))
  sigma.sq[r] <- ( 1/N)*( sum(term))
  rm( term) ;gc( verbose=FALSE)
  #message( "sigmasq of current micro iter", sigma.sq[r] )
  
  #level-2 variance component
  term <- list()
  term[[1]] <- tcrossprod( bi[[r]][[1]]) + 
    (	D[[r-1]] - 
        D[[r-1]] %*% t( Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
    )
  for(i in 2:n) 
    term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
    (	D[[r-1]] - 
        D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
    )
  term <- term[[n]]
  D[[r]] <- (1/n)*term
  rm( term) ;gc( verbose=FALSE)	 
  #message("D of current micro iter: ", D[[r]] )
  
  #level-1 and level-2 variance components (or typical or total variance)
  inv.Vi[[r]] <-list()
  for(i in 1:n){
    Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow = ni[i], ncol = ni[i])
    if(q==1)
      inv.Vi[[r]][[i]] <- 
        ( 1/sigma.sq[r]) * ( diag( rep(1,ni[i]))
                             -(( as.numeric(D[[r]])/sigma.sq[r])/( 1+ni[i]*( as.numeric(D[[r]])/sigma.sq[r])))*matrix( rep(1,(ni[i])^2)
                                                                                                                       , ncol = ni[i], nrow = ni[i]) )
    else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
  }
  Vi <- list(NULL) 
  inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
  
  
  #the generalized log-likelihood (GLL) 
  term <- vector( mode="numeric", length=n)
  for(i in 1:n)
    term[i] <- t( epsili[[i]]) %*% solve( sigma.sq[r]*diag( x = 1, nrow = ni[i], ncol = ni[i])) %*% epsili[[i]] +
    t( bi[[r]][[i]]) %*% solve( D[[r]]) %*% bi[[r]][[i]] +
    log( det(D[[r]])) +
    log( det(sigma.sq[r]*diag( x = 1, nrow = ni[i], ncol = ni[i])))
  GLL[r] <- sum( term)
  rm( term)
  gc( verbose=FALSE)
  
  #update the value of the Jump in GLL
  Jump[r] <- abs(( GLL[r]- GLL[r-1])/GLL[r] )
  
  if( Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
    convergence.iter[r] <- r
    if( verbose) message( "Converg. at iter no: ", r)
  } 
  
  
  ####################
  ####	STEP 3	####         
  ####################
  
  
  ###################################################
  #Iterating "F.niter" times to avoid early stopping#
  ###################################################
  F.niter=10
  I=1
  while(I<=F.niter){#repeat step 1 and 2
    
    ####################
    ####	STEP 1	####        
    ####################
    
    #update iteration number r
    r <- r+1
    if (verbose) 	
      message("MERF iter no: ", r)
    
    #update the length of the different lists
    
    length(sigma.sq) <- r
    length(D) <- r
    length(inv.Vi) <- r
    length(bi) <- r
    length(ci) <- r
    length(GLL) <- r
    
    length(Jump) <- r
    length(convergence.iter) <- r
    
    #update the transformed outcome, star.Yi
    star.Yi <- list() 
    for(i in 1:n){
      star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	- Xi[[i]] %*% ci[[r-1]]
    }
    
    #one STD random forest
    ######################
    
    for(i in 1:n)
      MERF.lDB[which(MERF.lDB$cluster.id==levels(MERF.lDB$cluster.id)[i]),]$star.Yi <- star.Yi[[i]] 
    rm(star.Yi) ;  gc(verbose=FALSE)
    
   
    fit.rf <- cforestmt( formula = fit.rf.formula , data = MERF.lDB	,weight_variable=weight_variable,  mtry = mtry, ntree = ntree, control = ctree_control(teststat="quad", testtype="Univ", mincriterion=0.95, minsplit = 50, minbucket = 50, maxsurrogate = min(1, ncol( xnam))), threads = threads)
    
    
    MERF.lDB$f.pred  <-  fit.rf$fitted$`(response)`
    
    #in matrix format
    fixed.pred <- list()	
    fixed.pred <- split(MERF.lDB, MERF.lDB$cluster.id) 	
    for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
    
    #fixed part
    ###########
    A = t(Xi[[1]]) %*% inv.Vi[[r-1]][[1]] %*% Xi[[1]]
    B = t(Xi[[1]]) %*% inv.Vi[[r-1]][[1]] %*% ( Yi[[1]] - fixed.pred[[1]])
    for(i in 2:n)
    { A= A + t(Xi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% Xi[[i]]
    B= B + t(Xi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% ( Yi[[i]] - fixed.pred[[i]])
    }
    ci[[r]] <- ginv(A) %*% B
    
    ci[r-1] <-list(NULL)
    #random	part
    ############
    for(i in 1:n)
      bi[[r]][[i]] <- D[[r-1]] %*% t( Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% ( Yi[[i]] - fixed.pred[[i]] - Xi[[i]] %*% ci[[r]])
    bi[r-1] <- list(NULL)		 
    
    
    ####################
    ####	STEP 2	####        
    ####################
    
    #level-1 variance component
    #residuals
    epsili <- list()
    for(i in 1:n)
      epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]] - Xi[[i]] %*% ci[[r]]
    
    
    term <- vector(mode="numeric",length=n)
    for(i in 1:n)
      term[i] <- crossprod(epsili[[i]]) + 
      sigma.sq[r-1] * (ni[i] - sigma.sq[r-1]* sum(diag(inv.Vi[[r-1]][[i]])))
    sigma.sq[r] <- (1/N)*(sum(term))
    rm(term) ;gc(verbose=FALSE)
    #message("sigmasq of current micro iter", sigma.sq[r] )
    
    #level-2 variance component
    term <- list()
    term[[1]] <- tcrossprod(bi[[r]][[1]]) + 
      (	D[[r-1]] - 
          D[[r-1]] %*% t(Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
      )
    for(i in 2:n) 
      term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
      (	D[[r-1]] - 
          D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
      )
    term <- term[[n]]
    D[[r]] <- (1/n)*term
    rm(term) ;gc(verbose=FALSE)	 
    #message("D of current micro iter: ", D[[r]] )
    
    #level-1 and level-2 variance components (or typical or total variance)
    inv.Vi[[r]] <-list()
    for(i in 1:n){
      Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
      if(q==1)
        inv.Vi[[r]][[i]] <- 
          (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                             -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                   , ncol=ni[i], nrow=ni[i]) )
      else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
    }
    Vi <- list(NULL) 
    inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
    
    
    #the generalized log-likelihood (GLL) 
    term <- vector(mode="numeric",length=n)
    for(i in 1:n)
      term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]] +
      t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]] +
      log(det(D[[r]])) +
      log(det(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
    GLL[r] <- sum(term)
    rm(term)
    gc(verbose=FALSE)
    
    #update the value of the Jump in GLL
    
    Jump[r] <- abs( (GLL[r]- GLL[r-1])/GLL[r] )
    
    if(Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
      convergence.iter[r] <- r
      if (verbose) message("Converg. at iter no: ", r)
      break
    } 
    I=I+1
    
  }
  #end for (I in 1: F.niter)
  ###########################
  
  ######################################
  #Iterating "max.niter" times at most #
  ######################################
  
  ###	END OF STEP 3	####        
  ############################
  
  
  #output to be returned (MERF model is the one at the last iteration)
  ###################### 
  
  
  output <- list(
    Jump[r]
    ,GLL
    ,convergence.iter
    ,fit.rf 
    ,bi[[r]]
    ,ci[[r]]
    ,sigma.sq[r]
    ,D	#,D[[r]]
  )
  
  names(output) <- c(
    "Jump[r]"
    ,"GLL"
    ,"convergence.iter"
    ,"fit.rf"
    ,"bi[[r]]"
    ,"ci[[r]]"
    ,"sigma.sq[r]"
    ,"D"	#,"D[[r]]"
  )
  
  #clean memory
  #############
  rm(
    xnam,
    MERF.lDB	
    ,ni,n,N
    ,Zi,q
    ,Yi	
    ,ntree
    ,mtry
    ,nodesize
    ,fit.rf.formula 
    ,fit.rf
    ,fixed.pred ,epsili
    
    ,sigma.sq,D,bi,Vi,inv.Vi
    
    ,F.niter ,max.niter
    ,smallest.Jump.allowed ,GLL ,Jump ,convergence.iter
    ,r,i,I
    
    ,verbose 		
  )
  gc(verbose=FALSE)
  
  
  #return
  #######
  output
  
}#end of MERF function

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



############################################################################################################
############################################################################################################
############################################################################################################
### end of the script	###
#sink()