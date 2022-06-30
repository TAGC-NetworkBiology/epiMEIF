library(randomForest)
library(dplyr)
library(varhandle)
Var_Names_Find <- function(t, index)
{
  ch <- sapply(index, function(i)unfactor(t$`split var`[i]))
  if(!any(is.na(ch)))
    return(ch)
}

trail_parents <- function(t, index)
{
  
 root=index
  while(root!= t$rowname[1])
  {
      if(index %in% t$`left daughter`)
        index=t$rowname[which(index==t$`left daughter`)]
      else
        index=t$rowname[which(index==t$`right daughter`)]
      root <- append(index, root)
  }
 return(root[-length(root)])
}  



ExtractInteraction<- function(i, t)
{
 
  parent_node <- vector()
  interaction_list <- list()
  if(t$status[i]!=-1)
  {
  parent_node <- trail_parents(t, i) 
  index_to_append <- list(c(as.numeric(t[i,1]), as.numeric(t[i,2])), c(as.numeric(t[i,1]), as.numeric(t[i,3])))
  interaction_list <- append(interaction_list, 
                             lapply(index_to_append, function(j)Var_Names_Find(t,j )))
  
  
  new_parent_node <- parent_node
  while(length(new_parent_node)!=0)
  { 
    index_to_append <- lapply(index_to_append, function(i) c(as.numeric(tail(new_parent_node,1)),i))
    interaction_list <- append(interaction_list, 
                               lapply(index_to_append, function(j)Var_Names_Find(t,j )))
    
    new_parent_node <- new_parent_node[-length(new_parent_node)]
  }
 
  return(interaction_list)  
  }
}

GetInteractionList <- function(rf, tree_index)
{
  t <- randomForest::getTree(rf, k = tree_index, labelVar = TRUE)  %>% tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
  if(any(t$status==-3))
  t$status[t$status==-3]=1
  interaction_list <- list()
  i=1
  if(i==1 & t$status==1)
  interaction_list <- append(interaction_list, ExtractInteraction(i, t))
  node_untraversed = vector()
  t$traverse=0
  while(i < nrow(t)&any(t$traverse[which(t$status==1)]!=2))
  {
      
      i_left=t$`left daughter`[i]
      i_right=t$`right daughter`[i]
      if(t$status[i_left]==-1 & t$status[i_right]==1)
        {t$traverse[i]=2
         i=i_right
         interaction_list <- append(interaction_list, ExtractInteraction(i, t))
        }
      else if(t$status[i_left]==1 & t$status[i_right]==-1) 
        {t$traverse[i]=2
         i=i_left
         interaction_list <- append(interaction_list, ExtractInteraction(i, t))
      }
      else if(t$status[i_left]==1 & t$status[i_right]==1) 
      {
        t$traverse[i]=1
        interaction_list <- append(interaction_list, ExtractInteraction(i_left, t))
        interaction_list <- append(interaction_list, ExtractInteraction(i_right, t))
        if((t$`right daughter`[i])<=nrow(t))
          node_untraversed =c(node_untraversed, t$`right daughter`[i])
        i=t$`left daughter`[i]
             }
      else
      {
        if(t$traverse[i]!=2)
        {interaction_list <- append(interaction_list, ExtractInteraction(i, t))
        t$traverse[i]=2
        }
        if(length(node_untraversed)!=0)
        {i=head(node_untraversed,1)
         t$traverse[match(i, t$`right daughter`)]=2
         node_untraversed = setdiff(node_untraversed, i)
        }  
      }
      
      
      
       #print(i)
      #print(node_untraversed)
    }
  interaction_list = sapply(interaction_list, function(ch) if(length(unique(ch))>1)return(unique(ch)))
  interaction_list = interaction_list[-which(sapply(interaction_list, is.null))]
  interaction_list =  unique(interaction_list)
  return(interaction_list)
  
  }
  






#row.matches <- function(y, X) {
#  i <- seq(nrow(X))
#  j <- 0
#  while(length(i) && (j <- j + 1) <= ncol(X)) 
#    i <- i[X[i, j] == y[j]]                                                                                                                                                                         i <- i[X[i, j] == y[j]]
#  i
#}

rowMatch <- function(A,B) {
  # Rows in A that match the rows in B
  # The row indexes correspond to A
  f <- function(...) paste(..., sep=":")
  if(!is.matrix(B)) B <- matrix(B, 1, length(B))
  a <- do.call("f", as.data.frame(A))
  b <- do.call("f", as.data.frame(B))
  match(b, a)
}
getInteractionMatrix_RF <- function(rf)
{
  
  ntree = rf$ntree
  Var_Interaction <- Var_Interaction_Mat <- list()
  for(ind in 1:ntree)
  {
    Var_Interaction[[ind]] <- GetInteractionList(rf, ind)
    Var_Interaction_Mat[[ind]] <- plyr::ldply(Var_Interaction[[ind]], rbind)
  } 
  
  All_Interactions <- do.call(qpcR:::rbind.na,  Var_Interaction_Mat)
  All_Interactions_Full <- Var_Interaction_Mat
  
  if(ncol(All_Interactions)==0)
    return(All_Interactions)
  else
  {
    colnames(All_Interactions) <- paste("Node", 1:ncol(All_Interactions), sep="")
    All_Interactions <- unique(All_Interactions)
    #All_Interactions <- All_Interactions%>% distinct()
    row_all_na <- sapply(1:nrow(All_Interactions), function(row)ifelse( all(is.na(All_Interactions[row,])), row,NA))
    if(length(na.omit(row_all_na))>0)
      All_Interactions <-All_Interactions[-na.omit(row_all_na),]
    
    
    
    All_Interactions[is.na(All_Interactions)]="Var_NA"
    All_Interactions_New <- t(apply(All_Interactions, 1, grr::sort2))
    All_Interactions_New <-  unique(All_Interactions_New)
    All_Interactions_Stats_New <- cbind( All_Interactions_New, matrix(0, ncol=ntree, nrow=nrow(All_Interactions_New), dimnames = list(paste("Interaction_",1:nrow(All_Interactions_New), sep=""),c( paste("Tree_", 1:ntree, sep="")))))
    
    
    system.time(for(i in  1:length(Var_Interaction))
    {
      if(nrow(All_Interactions_Full[[i]])!=0)
      {
        
        colnames(All_Interactions_Full[[i]])<-paste("V",1:ncol(All_Interactions_Full[[i]]), sep="_")
        All_Interactions_Full[[i]][is.na(All_Interactions_Full[[i]])]="Var_NA"
        All_Interactions_Full[[i]] <-  t(apply(All_Interactions_Full[[i]], 1, grr::sort2))
        #occ <- sapply(1:nrow(All_Interactions_New), function(interaction) sum(!is.na(row.match(All_Interactions_New[interaction,1:ncol(All_Interactions_Full[[i]])], All_Interactions_Full[[i]]))))
        occ <-  rowMatch(All_Interactions_Full[[i]], All_Interactions_New[,1:ncol(All_Interactions_Full[[i]])])
        occ[!is.na(occ)]=1
        occ[is.na(occ)]=0
        All_Interactions_Stats_New[,paste("Tree_", i, sep="")] =occ
      }
    })
    
    All_Interactions_Stats_New[All_Interactions_Stats_New=="Var_NA"]=""
    #All_Interactions_Stats[,paste("Tree_", 1:ntree, sep="")][ which( All_Interactions_Stats[, paste("Tree_", 1:ntree, sep="")]>=1, arr.ind=TRUE)]=1
    All_Interactions_Stats_New <- as.data.frame(All_Interactions_Stats_New)
    All_Interactions_Stats_New[, paste( "Tree_",1:ntree, sep="")]<- lapply( All_Interactions_Stats_New[, paste( "Tree_",1:ntree, sep="")], as.numeric)
    
    All_Interactions_Stats_New$Forest_Score <- rowSums( All_Interactions_Stats_New[, paste( "Tree_",1:ntree, sep="")])/ntree
    All_Interactions_Stats_New <- All_Interactions_Stats_New[ order( All_Interactions_Stats_New$Forest_Score),]
    All_Interactions_Stats_New <- All_Interactions_Stats_New[,-match(paste("Tree_",1:ntree,sep=""), colnames(All_Interactions_Stats_New))]
    
    return( All_Interactions_Stats_New)
  }
}

GenerateInteractionList_RF <- function(Interaction_List, Importance_Score_list, counter)
{
  
  All_Interactions_Stats <- do.call(qpcR:::rbind.na, Interaction_List)
  All_Interactions_Stats <-All_Interactions_Stats %>% distinct()
  All_Interactions_Stats <- do.call(qpcR:::rbind.na, sapply(1:nrow(All_Interactions_Stats), function(row)All_Interactions_Stats[row,order(as.numeric(do.call(rbind, strsplit(unlist(All_Interactions_Stats[row,]), "_"))[,2]))]))
  colnames(All_Interactions_Stats) <- paste("Node",1:ncol(All_Interactions_Stats), sep="")
  All_Interactions_Stats <-All_Interactions_Stats %>% distinct()
  
  #for(j in 1: nrow(All_Interactions_Stats_SRF))
  #{
  #  print(j)
  #  str <- All_Interactions_Stats_SRF[j,]
  #  str[which(str=="")]=NA
  #  str <-  str[!is.na(str)]
  #  for(k in (j+1):nrow(All_Interactions_Stats_SRF))
  #  {
  #    strnew <- All_Interactions_Stats_SRF[k,]
  #    strnew[which(strnew=="")]=NA
  #    strnew <-  strnew[!is.na(strnew)]
  #    if(length(str)==length(strnew)&all(str %in% strnew))
  #      {All_Interactions_Stats_SRF <- All_Interactions_Stats_SRF[-k,]
  #       k=k-1
  #      }
  #  }
  #}
  
  All_Interactions <- All_Interactions_Stats
  All_Interactions_Stats <- cbind(All_Interactions_Stats, matrix(0, ncol=NiterChild, nrow=nrow(All_Interactions_Stats), dimnames = list(paste("Interaction_",1:nrow(All_Interactions_Stats), sep=""),c( paste("Forest_", 1:NiterChild, sep="")))))
  All_Interactions[All_Interactions==""]=NA
  system.time(for(i in 1:NiterChild)
  { if(counter=="SRF")
    mat_to_check <- Importance_Score_list[[i]][[3]]

  mat_to_check <- mat_to_check%>% filter(Forest_Score>0.04)
  mat_to_check <- cbind(do.call(qpcR:::rbind.na, sapply(1:nrow(mat_to_check), function(row)mat_to_check[,-ncol(mat_to_check)][row,order(as.numeric(do.call(rbind, strsplit(unlist(mat_to_check[,-ncol(mat_to_check)][row,]), "_"))[,2]))])), Forest_Score=mat_to_check$Forest_Score)
  
  #mat_to_check[mat_to_check ==""]=NA
  #pos_match <- sapply(1:nrow(All_Interactions_Stats_SRF), function(row) ifelse(is.na(All_Interactions_Stats[row,setdiff(colnames(mat_to_check), "Forest_Score")], mat_to_check[,-ncol(mat_to_check)])), 0, mat_to_check$Forest_Score[row.match(All_Interactions_Stats[row,setdiff(colnames(mat_to_check), "Forest_Score")], mat_to_check[,-ncol(mat_to_check)])])) 
  #pos_match <-sapply(1:nrow(All_Interactions_Stats),function(row) na.omit(RowContain(All_Interactions[row,], mat_to_check[,-ncol(mat_to_check)])))
  #pos_match <- sapply(pos_match, function(pos) ifelse(length(pos)==0, NA, pos))
  pos_match_new <-sapply(1:nrow(All_Interactions_Stats),function(row) na.omit(row.match(All_Interactions[row,1:(ncol(mat_to_check)-1)], mat_to_check[,-ncol(mat_to_check)])))
  pos_match_new <- sapply(pos_match_new, function(pos) ifelse(length(pos)==0, NA, pos))
  All_Interactions_Stats[[paste("Forest_",i,sep="")]] <- mat_to_check$Forest_Score[pos_match_new]
  #All_Interactions_Stats[[paste("Forest_",i,sep="")]][All_Interactions_Stats[[paste("Forest_",i,sep="")]]>1]=1
  })
  
  All_Interactions_Stats[All_Interactions_Stats==""]=NA
  All_Interactions_Stats$Sum_Forest_Score <- sapply( 1:nrow(All_Interactions_Stats), function(i) sum( as.numeric(All_Interactions_Stats[ i, paste("Forest_", 1:10,sep="")]), na.rm=TRUE) )
  All_Interactions_Stats$Median_Forest_Score <- sapply( 1:nrow(All_Interactions_Stats), function(i)median( as.numeric(All_Interactions_Stats[ i, paste("Forest_", 1:10,sep="")]), na.rm=TRUE) )
  
  All_Interactions_Stats<- All_Interactions_Stats%>% distinct()
  All_Interactions_Stats<- All_Interactions_Stats[order(All_Interactions_Stats$Median_Forest_Score),]
  All_Interactions_Stats<-All_Interactions_Stats[which(All_Interactions_Stats$Median_Forest_Score>=0.08),]
  All_Interactions_Stats$Rank <- rank(max(All_Interactions_Stats$Median_Forest_Score)-All_Interactions_Stats$Median_Forest_Score)
  
  return(All_Interactions_Stats)
}
  
GenerateInteractions_RF <- function(Interaction_List, Importance_Score_list, counter)
{
  
  All_Interactions_Stats <- do.call(qpcR:::rbind.na, Interaction_List)
  All_Interactions_Stats <-All_Interactions_Stats %>% distinct()
  All_Interactions_Stats <- do.call(qpcR:::rbind.na, sapply(1:nrow(All_Interactions_Stats), function(row)All_Interactions_Stats[row,order(as.numeric(do.call(rbind, strsplit(unlist(All_Interactions_Stats[row,]), "_"))[,2]))]))
  colnames(All_Interactions_Stats) <- paste("Node",1:ncol(All_Interactions_Stats), sep="")
  All_Interactions_Stats <-All_Interactions_Stats %>% distinct()
  
  #for(j in 1: nrow(All_Interactions_Stats_SRF))
  #{
  #  print(j)
  #  str <- All_Interactions_Stats_SRF[j,]
  #  str[which(str=="")]=NA
  #  str <-  str[!is.na(str)]
  #  for(k in (j+1):nrow(All_Interactions_Stats_SRF))
  #  {
  #    strnew <- All_Interactions_Stats_SRF[k,]
  #    strnew[which(strnew=="")]=NA
  #    strnew <-  strnew[!is.na(strnew)]
  #    if(length(str)==length(strnew)&all(str %in% strnew))
  #      {All_Interactions_Stats_SRF <- All_Interactions_Stats_SRF[-k,]
  #       k=k-1
  #      }
  #  }
  #}
  
  All_Interactions <- All_Interactions_Stats
  All_Interactions_Stats <- cbind(All_Interactions_Stats, matrix(0, ncol=NiterChild, nrow=nrow(All_Interactions_Stats), dimnames = list(paste("Interaction_",1:nrow(All_Interactions_Stats), sep=""),c( paste("Forest_", 1:NiterChild, sep="")))))
  All_Interactions[All_Interactions==""]=NA
  system.time(for(i in 1:NiterChild)
  { if(counter=="SRF")
    mat_to_check <- Importance_Score_list[[i]]
  
  mat_to_check <- mat_to_check%>% filter(Forest_Score>0.04)
  mat_to_check <- cbind(do.call(qpcR:::rbind.na, sapply(1:nrow(mat_to_check), function(row)mat_to_check[,-ncol(mat_to_check)][row,order(as.numeric(do.call(rbind, strsplit(unlist(mat_to_check[,-ncol(mat_to_check)][row,]), "_"))[,2]))])), Forest_Score=mat_to_check$Forest_Score)
  
  #mat_to_check[mat_to_check ==""]=NA
  #pos_match <- sapply(1:nrow(All_Interactions_Stats_SRF), function(row) ifelse(is.na(All_Interactions_Stats[row,setdiff(colnames(mat_to_check), "Forest_Score")], mat_to_check[,-ncol(mat_to_check)])), 0, mat_to_check$Forest_Score[row.match(All_Interactions_Stats[row,setdiff(colnames(mat_to_check), "Forest_Score")], mat_to_check[,-ncol(mat_to_check)])])) 
  #pos_match <-sapply(1:nrow(All_Interactions_Stats),function(row) na.omit(RowContain(All_Interactions[row,], mat_to_check[,-ncol(mat_to_check)])))
  #pos_match <- sapply(pos_match, function(pos) ifelse(length(pos)==0, NA, pos))
  pos_match_new <-sapply(1:nrow(All_Interactions_Stats),function(row) na.omit(row.match(All_Interactions[row,1:(ncol(mat_to_check)-1)], mat_to_check[,-ncol(mat_to_check)])))
  pos_match_new <- sapply(pos_match_new, function(pos) ifelse(length(pos)==0, NA, pos))
  All_Interactions_Stats[[paste("Forest_",i,sep="")]] <- mat_to_check$Forest_Score[pos_match_new]
  #All_Interactions_Stats[[paste("Forest_",i,sep="")]][All_Interactions_Stats[[paste("Forest_",i,sep="")]]>1]=1
  })
  
  All_Interactions_Stats[All_Interactions_Stats==""]=NA
  All_Interactions_Stats$Sum_Forest_Score <- sapply( 1:nrow(All_Interactions_Stats), function(i) sum( as.numeric(All_Interactions_Stats[ i, paste("Forest_", 1:10,sep="")]), na.rm=TRUE) )
  All_Interactions_Stats$Median_Forest_Score <- sapply( 1:nrow(All_Interactions_Stats), function(i)median( as.numeric(All_Interactions_Stats[ i, paste("Forest_", 1:10,sep="")]), na.rm=TRUE) )
  
  All_Interactions_Stats<- All_Interactions_Stats%>% distinct()
  All_Interactions_Stats<- All_Interactions_Stats[order(All_Interactions_Stats$Median_Forest_Score),]
  All_Interactions_Stats<-All_Interactions_Stats[which(All_Interactions_Stats$Median_Forest_Score>=0.08),]
  All_Interactions_Stats$Rank <- rank(max(All_Interactions_Stats$Median_Forest_Score)-All_Interactions_Stats$Median_Forest_Score)
  
  return(All_Interactions_Stats)
}



