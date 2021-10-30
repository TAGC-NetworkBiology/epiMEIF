WORKING_DIR <-  "D:/workspace/MERF/output_files"
INPUT_DIR <- NEW_DIR <- "D:/workspace/From_Meso"
DATA_FILE <- file.path(WORKING_DIR, "Heartperiod_Median_Data.RDS")
#Phenotype_Data <- readRDS( DATA_FILE)
#saveRDS(Genotype_Data,"Data_Genotype_NMiss.RDS")
#Genotype_Data <- readRDS(file.path(NEW_DIR, "Data_Genotype_Imputed.RDS"))
#dim(Genotype_Data)

Phenotype_Data <- readRDS( DATA_FILE)
Genotype_Data <- readRDS(file.path(INPUT_DIR, "Data_Genotype_Imputed.RDS"))
ALT_SCRIPT_DIR <- "D:/workspace/MERF/output_files/RUN_2020_11_06/script/03_post_gwas_analysis"
source(file.path(ALT_SCRIPT_DIR, "Interaction_Score_partykit-v4.R"))
source(file.path(ALT_SCRIPT_DIR, "Simulation_Functions_Age1-v2.R"))

###Extracting the first 1000 SNPS
Genotype_Data_Sub <- Genotype_Data[,1:1000]
rownames(Genotype_Data_Sub) <- gsub("line_","dgrp",rownames(Genotype_Data))
strain <-rownames(Genotype_Data_Sub)
Genotype_Data_Sub <- as_tibble(Genotype_Data_Sub)
Genotype_Data_Sub$strain <- strain


###Computing the MAF of each SNP
MAF <- sapply(1:(ncol(Genotype_Data_Sub)-1), function(i) min(table(Genotype_Data_Sub[,i])/sum(table(Genotype_Data_Sub[,i]))))

###Simulating the SNPs based on the MAF of each SNP
set.seed(1000)
Genotype_Data_Simulated <- sapply(1:(ncol(Genotype_Data_Sub)-1), function(i)rbern(167,1-MAF[i]))
colnames(Genotype_Data_Simulated) <- names(MAF) <- paste("SNP_",1:1000,sep="")

Genotype_Data_Simulated[Genotype_Data_Simulated==0]=2
Genotype_Data_Simulated[Genotype_Data_Simulated==1]=0

Genotype_Data_Simulated <- cbind(Genotype_Data_Simulated,strain= Genotype_Data_Sub$strain)

Genotype_Data_Simulated[,setdiff( colnames(Genotype_Data_Simulated), "strain")] <- sapply(Genotype_Data_Simulated[,setdiff( colnames(Genotype_Data_Simulated), "strain")], as.factor)

Pheno_Geno_Data <- merge(Phenotype_Data, as.data.frame(Genotype_Data_Simulated), by="strain", all.x=TRUE)

Pheno_Geno_Data_Age1 <- Pheno_Geno_Data %>% filter(age==1)
Pheno_Geno_Data_Age1[colnames(Pheno_Geno_Data_Age1)[8:1007]] <- lapply(Pheno_Geno_Data_Age1[colnames(Pheno_Geno_Data_Age1)[8:1007]] , as.factor)
Pheno_Geno_Data_Age1_Orig <- Pheno_Geno_Data_Age1




MAF <- MAF[which(MAF > 0.1&MAF!=1)]
Genotype_Data_Simulated <- Genotype_Data_Simulated[ , match(names(MAF), colnames(Genotype_Data_Simulated))]
Genotype_Data_Simulated <- as.data.frame(Genotype_Data_Simulated)
Genotype_Data_Simulated[setdiff( colnames(Genotype_Data_Simulated), "strain")] <- sapply(Genotype_Data_Simulated[setdiff( colnames(Genotype_Data_Simulated), "strain")], as.numeric)

set.seed(100)
Seed_Root <- readRDS("D:/workspace/From_Meso/Simulation_Scenario_Codes_Meso/Final_Results/Seed_Aging_Scenario1_30interactors.RDS")
NiterRoot <- length(Seed_Root)
NiterChild=10


combn_mat <- combn(1:8,2)
SNP_Pairs_Significance <- matrix( NA, nrow=length(Seed_Root), ncol=ncol(combn_mat), dimnames=list(Seed_Root
                                                                                                  , paste("Interaction_",1:ncol(combn_mat),sep="")))
SNP_Significance <- matrix( NA, nrow=length(Seed_Root), ncol=8, dimnames=list(Seed_Root, paste("SNP_",1:8,sep="")))



TopInteractionsCount_v2<-  matrix(0, nrow=NiterRoot, ncol=ncol(combn_mat)*2+2, dimnames = list(1:NiterRoot, c(sapply(1:ncol(combn_mat),function(col)paste("Simble_Interaction", paste(combn_mat[,col], collapse = ":"), sep="")),sapply(1:ncol(combn_mat),function(col)paste("Weighted_Interaction", paste(combn_mat[,col], collapse = ":"), sep="")), "Interaction_validated_SRF", "Interaction_validated_WRF")))


seed_root=Seed_Root[1]

  Pheno_Geno_Data_Age1=Pheno_Geno_Data_Age1_Orig
  index=match(seed_root, Seed_Root)
  print(match(seed_root, Seed_Root))
  set.seed(seed_root)
  SNPs_Selected <- sample(names(MAF)[1:100],8)
  
  
  #Genotype_Corr <- cor(Genotype_Data_Simulated[,match(SNPs_Selected, colnames(Genotype_Data_Simulated))], use="pairwise.complete.obs", method="spearman")
  Pheno_Geno_Data_Age1_simulation <- Pheno_Geno_Data_Age1[, c(colnames(Pheno_Geno_Data_Age1)[1:7], SNPs_Selected)]
  sigma_beta <- sqrt(0.375)
  sigma_v <- sqrt(0.125)
  set.seed(1000)
  beta<- rnorm(8,0, sigma_beta)
  set.seed(1000)
  err <- rnorm(nrow(Pheno_Geno_Data_Age1_simulation), 0, sigma_v)
  Pheno_Geno_Data_Age1_simulation[SNPs_Selected] <- sapply(Pheno_Geno_Data_Age1_simulation[SNPs_Selected], as.numeric)-1
  
  Pheno_Geno_Data_Age1_simulation$Y_new <-rowSums( c(beta[1], beta[2], beta[3], beta[4])*Pheno_Geno_Data_Age1_simulation[,SNPs_Selected[1:4]])+
    beta[7]*0.7*Pheno_Geno_Data_Age1_simulation[,SNPs_Selected[5]]*Pheno_Geno_Data_Age1_simulation[,SNPs_Selected[6]]+
    +beta[8]*Pheno_Geno_Data_Age1_simulation[,SNPs_Selected[7]]*Pheno_Geno_Data_Age1_simulation[,SNPs_Selected[8]]+err
  
  Pheno_Geno_Data_Age1$Y_new <- Pheno_Geno_Data_Age1_simulation$Y_new
  
  SingleMarkerSignificance <- unlist(sapply( colnames(Genotype_Data_Simulated), function(snp){if(length(table(Genotype_Data_Simulated[,snp]))>1){t <- anova(lm(as.formula( paste("Y_new~", snp, sep=""))  , data=Pheno_Geno_Data_Age1)); return( t$`Pr(>F)`[1])}}))
  SingleMarkerSignificance_SNPs_Selected <- SingleMarkerSignificance[match(SNPs_Selected, names(SingleMarkerSignificance))]
  SNP_Significance[index, ] <- SingleMarkerSignificance_SNPs_Selected
  SNPs_for_cf <- setdiff(names(which(SingleMarkerSignificance>0.1)), SNPs_Selected)[1:50]
  Pheno_Geno_Data_Age1 <- Pheno_Geno_Data_Age1[, c("Y_new", SNPs_Selected, SNPs_for_cf)] 
  Pheno_Geno_Data_Age1[c(SNPs_Selected, SNPs_for_cf)]<- lapply(Pheno_Geno_Data_Age1[c(SNPs_Selected, SNPs_for_cf)], factor)
  
  #fit_removed <- lm(as.formula(paste("Y_new~", paste(sapply(setdiff(1:ncol(combn_mat), c(23,28)), 
  #                                           function(row)paste(SNPs_Selected[combn_mat[,row]], collapse=":")), 
  #                                 collapse="+"), sep="")), data=Pheno_Geno_Data_Age1_simulation)
  #Pheno_Geno_Data_Age1_simulation$Y_new2 <- Pheno_Geno_Data_Age1_simulation$Y_new-predict(fit_removed)
  #par(mfrow=c(1,2))
  #boxplot(Pheno_Geno_Data_Age1_simulation$Y_new)  
  #boxplot(Pheno_Geno_Data_Age1_simulation$Y_new2)  
  #Pheno_Geno_Data_Age1$Y_new <- Pheno_Geno_Data_Age1_simulation$Y_new
  colnames(Pheno_Geno_Data_Age1)[1]="PHENOTYPE"
  #Pheno_Geno_Data_Age1 <- Pheno_Geno_Data_Age1[, -c(1:7)]
  saveRDS(Pheno_Geno_Data_Age1, "D:/workspace/From_Meso/SingleGWASFiles/GWAS_Final/R_Codes_Final/Age1_Dataset.RDS")
  write.csv(Pheno_Geno_Data_Age1, "D:/workspace/From_Meso/SingleGWASFiles/GWAS_Final/R_Codes_Final/Age1_Dataset.csv")  
   