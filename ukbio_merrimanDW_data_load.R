library(RMySQL)
library(data.table)
library(dplyr)
drv = dbDriver("MySQL")
db = dbConnect(drv, default.file = '~/.my.cnf', dbname="merrimanDW")

dbGetQuery(db,'show tables;')

dbGetQuery(db, 'describe dimSourceDataset;')
dbGetQuery(db, 'select * from dimVar')
dbGetQuery(db, 'select * from factGWAS')
dbGetQuery(db, 'select * from dimExperiment')
dbGetQuery(db, 'describe dimExperiment;')
dbGetQuery(db, 'describe dimCondition;')

scratch_dir <- "/media/xsan/scratch/merrimanlab/murray/working_dir/UkBio/"

#ukbio dataset has been already loaded
datasetID <- dbGetQuery(db, 'select datasetID FROM dimSourceDataset WHERE datasetName = "UKBIOBANK";')



# dimCondition
conditions <- c('all', 'hospital', 'winnard','self report', 'self report and ULT', 'all (male)', 'hospital (male)', "no diuretic", "diuretic")
for(cond in conditions){
  dbGetQuery(db, paste0('insert into dimCondition (conditionDescription) values', '("',cond,'")'))
}


g_cond <- c('all', 'hospital', 'winnard','self report', 'self report and ULT', 'all (male)', 'hospital (male)')
c_cond <- c('all',"no diuretic", "diuretic")

control_switch <- function(type){
  switch(type, 
         all = 'controls',
         `no diuretic` = 'controls_no_diuretics',
         diuretic = 'controls_diuretics'
  )
}
case_switch <-function(type){
  switch(type,
         all = 'all',
         winnard = 'winnard',
         `self report` = 'self',
         `self report and ULT` = 'self_ult',
         `all (male)` = 'all_male',
         `hospital (male)` = 'hosp_male',
         hospital = 'hosp'
  )
}

add_dimExperiment <- function(gc,cc){
  # dimExperiment
  dbGetQuery(db, 'describe dimExperiment;')
  dimExperiment <- data.frame( trait = "Gout", 
                               datasetID = datasetID, 
                               Description = paste("Gout unadjusted UK Biobank using",gc,"cases",cc,"controls"), 
                               caseCondition = case_condID[,1], 
                               controlCondition = control_condID[,1], 
                               Covariates = "", 
                               ParameterType = 0) # 0 = OR, 1 = Beta
  
  
  dimExperimentCovar <- data.frame( trait = "Gout", 
                                    datasetID = datasetID, 
                                    Description = paste("Gout adjusted UK Biobank using",gc,"cases",cc,"controls"), 
                                    caseCondition = case_condID[,1], 
                                    controlCondition = control_condID[,1], 
                                    Covariates = "sex age waist waist_to_height_ratio", 
                                    ParameterType = 0) # 0 = OR, 1 = Beta
  
  dbWriteTable(db, name='dimExperiment',dimExperiment, append = TRUE,row.names = FALSE)
  dbWriteTable(db, name='dimExperiment',dimExperimentCovar, append = TRUE,row.names = FALSE)
}

fam_file_load<- function(gc,cc){
  message('read fam')
  fam_file <- fread(paste0(scratch_dir,'GWAS_all_controls/',control_switch(cc),'/',case_switch(gc),'/chrallimpv1.fam_',case_switch(gc)), header=FALSE)
  colnames(fam_file) <- c('fid','iid','PID','MID','SEX','AFF')
  
  fam_file[AFF == 1, .(fid, iid)]
  fam_file[AFF == 2, .(fid,iid)]
  fam_file$CaseControl <- fam_file$AFF
  return(fam_file)
}


factGWAS_dimVariant <- function(gc,cc,ex){
  #fills factGWAS and dimVariant
  exCovar <- dbGetQuery(db,paste0('select Covariates from dimExperiment where experimentID = ',ex ,';'))
  if(exCovar == ""){
    adj = 'unadjusted'
  } else {
    adj = 'adjusted'
  }
  message(adj)
  for( chr in 1:22){
    message(chr)
    message('read assoc')
    results <- fread(paste0(scratch_dir,'GWAS_all_controls/',control_switch(cc),'/',case_switch(gc),'/',adj,'/',control_switch(cc),case_switch(gc),'_chr',chr,'.assoc.logistic.tab'), header=TRUE, colClasses = c('integer','character','integer','character','character','integer',rep('numeric',6)))
    colnames(results)
    results$row <- as.numeric(rownames(results))
    results$obs <- results$NMISS
    results <- results[TEST == 'ADD' & SNP != '.']
    results$row <- as.numeric(rownames(results))
    setkey(results, SNP, obs)
    
    
    # row number ~~should~~ (observed but not proved) match for joining freq with the snp per row version of hwe
    message('read freq')
    freq <- fread(paste0(scratch_dir,'GWAS_all_controls/',control_switch(cc),'/',case_switch(gc),'/',adj,'/',control_switch(cc),case_switch(gc),'_chr',chr,'.frq.cc.tab'), header=TRUE)
    freq$row <- as.numeric(rownames(freq))
    freq$obs <- apply(freq[,.(NCHROBS_A,NCHROBS_U)],1,sum) / 2
    setkey(freq, row)
    colnames(freq) <- c("CHR_frq","SNP_frq","A1_frq","A2_frq","MAF_A","MAF_U","NCHROBS_A","NCHROBS_U","row","obs")      
    
    
    #use row to separate and rejoin hwe back together as single row per snp
    message('read hwe')
    hwe <- fread(paste0(scratch_dir,'GWAS_all_controls/',control_switch(cc),'/',case_switch(gc),'/',adj,'/',control_switch(cc),case_switch(gc),'_chr',chr,'.hwe.tab'), header=TRUE)
    hwe$obs_hwe <- unlist(lapply(lapply(strsplit(hwe$GENO, '/'), as.numeric),sum))
    hwe$row1 <- as.numeric(rownames(hwe))
    
    hwe_all <- hwe[TEST == 'ALL']
    hwe_all$row <- as.numeric(rownames(hwe_all))
    setkey(hwe_all, row)
    names(hwe_all) <- c("CHR_ALL", "SNP_ALL", "TEST_ALL", "A1_ALL", "A2_ALL", "HWE_GENO_ALL", "O(HET)_ALL", "E(HET)_ALL", "HWE_P_ALL", "obs_ALL","row1_all","row")
    
    hwe_a <- hwe[TEST == 'AFF']
    hwe_a$row <- as.numeric(rownames(hwe_a))
    setkey(hwe_a, row)
    names(hwe_a) <- c("CHR_A", "SNP_A", "TEST_A", "A1_A", "A2_A", "HWE_GENO_A", "O(HET)_A", "E(HET)_A", "HWE_P_A", "obs_A","row1_A","row")   
    
    hwe_u <- hwe[TEST == 'UNAFF']
    hwe_u$row <- as.numeric(rownames(hwe_u))
    setkey(hwe_u, row)
    names(hwe_u) <- c("CHR_U", "SNP_U", "TEST_U", "A1_U", "A2_U", "HWE_GENO_U", "O(HET)_U", "E(HET)_U", "HWE_P_U", "obs_U","row_U","row")  
    hwe_all_a <- hwe_all[hwe_a]
    hwe_au <- hwe_all_a[hwe_u]
    freq_hwe_au <- freq[hwe_au]
    freq_hwe_au  %>% filter(obs != obs_A + obs_U)
    freq_hwe_au$SNP <- freq_hwe_au$SNP_frq
    setkey(freq_hwe_au,SNP, obs)
    results_freq_hwe_au <- results[freq_hwe_au]
    results_freq_hwe_au  %>% filter(SNP != SNP_frq | SNP != SNP_A | SNP != SNP_U )
    results_freq_hwe_au  %>%filter(A1 != A1_frq | A1 != A1_A | A1 != A1_U)
    
    results_freq_hwe_au <- results_freq_hwe_au[!is.na(results_freq_hwe_au$BP)]
    results_freq_hwe_au <- results_freq_hwe_au[!(SNP %in% results_freq_hwe_au$SNP[duplicated(results_freq_hwe_au$SNP)])]
    results_freq_hwe_au$MAF <- unlist(lapply(lapply(strsplit(results_freq_hwe_au$HWE_GENO_ALL, '/'), as.numeric), function(x){(2* x[1] + x[2])/(2*sum(x))}))
    setkey(results_freq_hwe_au,row)
    rm(freq, hwe, hwe_all, results, hwe_a, hwe_u)
    
    
    #remove variants that are not SNPs
    dim(results_freq_hwe_au[sapply(results_freq_hwe_au$A2_A,nchar) == 1 | sapply(results_freq_hwe_au$A1,nchar) ==1 ])
    results_freq_hwe_au <- results_freq_hwe_au[sapply(results_freq_hwe_au$A2_A,nchar) == 1 & sapply(results_freq_hwe_au$A1,nchar) ==1 ]
    variants <- results_freq_hwe_au[, .(CHR, SNP, BP, A1, A2_frq, row)]
    colnames(variants) <- c("chromosome","snp","GRCh37_bp","A1","A2", 'row')
    setkey(variants, chromosome,snp,A1,A2)
    variants$GRCh38_bp <- NA
    dimVar <-  as.data.table(dbGetQuery(db, paste0('SELECT * FROM dimVariant WHERE chromosome = ',chr)))
    setkey(dimVar, chromosome,snp,A1,A2)
    dimVar_orig <- dimVar
    head(dimVar[variants])
    
    # add missing variants
    varToAdd <- dimVar[variants]
    if(nrow(varToAdd[(!is.na(GRCh37_bp) & GRCh37_bp != i.GRCh37_bp) | (!is.na(GRCh38_bp) & GRCh38_bp != i.GRCh38_bp)]) == 0){
      varToAdd <- varToAdd[is.na(varToAdd$variantID)]
      varToAdd$GRCh37_bp <- varToAdd$i.GRCh37_bp
      varToAdd$GRCh38_bp <- varToAdd$i.GRCh38_bp
      dbWriteTable(db, name='dimVariant',varToAdd[,-c('row',"i.GRCh37_bp","i.GRCh38_bp"),with=FALSE], append = TRUE,row.names = FALSE)
    }else{
      message("Mismatch in known positions for varToAdd")
    }
    
    #refresh dimVar to get all variantIDs
    dimVar <-  as.data.table(dbGetQuery(db, paste0('SELECT * FROM dimVariant WHERE chromosome = ',chr)))
    setkey(dimVar, chromosome,snp,A1,A2)
    
    if(sum(is.na(dimVar$variantID) > 0)){
      message("missing variantIDs")
    }
    
    
    #join variantIDs onto results
    
    variants2 <- variants[dimVar]
    setkey(variants2, row)
    
    variants3 <- variants2[!duplicated(variants2$row)]
    
    results_freq_hwe_au2 <- results_freq_hwe_au[variants3]
    
    
    factTable <- results_freq_hwe_au2[,.(SNP,OR,SE, P, MAF_A, MAF_U, MAF, HWE_GENO_A, HWE_GENO_U,HWE_GENO_ALL, HWE_P_A, HWE_P_U, HWE_P_ALL,variantID)]
    
    
    colnames(factTable) <- c("SNP","parameterEstimate","stdError","pValue","MAF_A","MAF_U","MAF_ALL","HWE_GENO_A", "HWE_GENO_U", "HWE_GENO_ALL","HWE_P_A", "HWE_P_U","HWE_P_ALL", 'variantID')
    factTable$majorAllele  <- 2
    factTable$experimentID <- ex
    
    
    dbWriteTable(db, name='factGWAS',na.omit(factTable[, -'SNP', with = FALSE]), append = TRUE,row.names = FALSE)
  }
}
system.time(
for(gc in g_cond){
  for(cc in c_cond){
    
    control_cond <- cc
    control_condID <- dbGetQuery(db, paste0('select conditionID from dimCondition where conditionDescription = "',control_cond,'";'))
    
    case_cond <- gc
    case_condID <- dbGetQuery(db, paste0('select conditionID from dimCondition where conditionDescription = "',case_cond,'";'))
    message((paste(gc,cc)))
    
    add_dimExperiment(gc = gc, cc = cc)
    
    expID <- dbGetQuery(db, paste0('select experimentID from dimExperiment WHERE caseCondition = ', '"',case_condID, '"',' and controlCondition = ', '"',control_condID, '"' ,' and datasetID =', datasetID,';'))
    #dimPerson
    fam_file <- fam_file_load(gc = gc, cc = cc)
    for(ex in expID$experimentID){
      message(paste(cc,gc,ex))
      fam_file$experimentID <-  ex
      dbWriteTable(db, name = 'dimPeople', fam_file[, .(fid,iid, experimentID, CaseControl)], append=TRUE,row.names = FALSE)
      #factGwas and dimVariant
      factGWAS_dimVariant(gc = gc,cc = cc,ex = ex)
    }
    
  }
}
)


