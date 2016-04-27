library(RMySQL)
library(data.table)
library(dplyr)
drv = dbDriver("MySQL")
db = dbConnect(drv, default.file = '~/.my.cnf', dbname="merrimanDW")

dbGetQuery(db,'show tables;')

dbGetQuery(db, 'describe dimSourceDataset;')
dbGetQuery(db, 'select * from dimExperiment')
dbGetQuery(db, 'describe dimExperiment;')
dbGetQuery(db, 'describe dimCondition;')

scratch_dir <- "/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/"

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

for(gc in g_cond){
  for(cc in c_cond){
    
    control_cond <- cc
    control_condID <- dbGetQuery(db, paste0('select conditionID from dimCondition where conditionDescription = "',control_cond,'";'))
    
    case_cond <- gc
    case_condID <- dbGetQuery(db, paste0('select conditionID from dimCondition where conditionDescription = "',case_cond,'";'))
    
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
    
    expID <- dbGetQuery(db, paste0('select experimentID from dimExperiment WHERE caseCondition = ', '"',case_condID, '"',' and controlCondition = ', '"',control_condID, '"' ,' and datasetID =', datasetID,';'))
    
       
    
    #dimPerson
    fam_file <- fread(paste0(scratch_dir,'GWAS_all_controls/',control_switch(cc),'/',case_switch(gc),'/chrallimpv1.fam_',case_switch(gc)), header=FALSE)
    colnames(fam_file) <- c('fid','iid','PID','MID','SEX','AFF')
    
    fam_file[AFF == 1, .(fid, iid)]
    fam_file[AFF == 2, .(fid,iid)]
    fam_file$CaseControl <- fam_file$AFF
    for(ex in expID$experimentID){
    fam_file$experimentID <-  ex
        dbWriteTable(db, name = 'dimPeople', fam_file[, .(fid,iid, experimentID, CaseControl)], append=TRUE,row.names = FALSE)
    }
    
  }
}


for( chr in 1:1){
results <- fread(paste0(scratch_dir,'GWAS_all_controls/',control_switch(cc),'/',case_switch(gc),'/adjusted/',control_switch(cc),case_switch(gc),'_chr',chr,'.assoc.logistic'), header=TRUE)
colnames(results)
results$row <- as.numeric(rownames(results))
setkey(results, row)
results <- results[TEST == 'ADD' & SNP != '.']



freq <- fread(paste0(scratch_dir,'GWAS_all_controls/',control_switch(cc),'/',case_switch(gc),'/adjusted/',control_switch(cc),case_switch(gc),'_chr',chr,'.frq.cc'), header=TRUE)
freq$row <- as.numeric(rownames(freq))
setkey(freq, row)
freq <- freq[SNP != '.']
colnames(freq) <- c("CHR_frq","SNP_frq","A1_frq","A2_frq","MAF_A","MAF_U","NCHROBS_A","NCHROBS_U","row")      

hwe <- fread(paste0(scratch_dir,'GWAS_all_controls/',control_switch(cc),'/',case_switch(gc),'/adjusted/',control_switch(cc),case_switch(gc),'_chr',chr,'.hwe'), header=TRUE)
hwe_all <- hwe[TEST == 'ALL' & SNP != '.']
hwe_all$row <- as.numeric(rownames(hwe_all))
setkey(hwe_all, row)
names(hwe_all) <- c("CHR_ALL", "SNP_ALL", "TEST_ALL", "A1_ALL", "A2_ALL", "HWE_GENO_ALL", "O(HET)_ALL", "E(HET)_ALL", "HWE_P_ALL", "row")

hwe_a <- hwe[TEST == 'AFF']
hwe_a$row <- as.numeric(rownames(hwe_a))
setkey(hwe_a, row)
names(hwe_a) <- c("CHR_A", "SNP_A", "TEST_A", "A1_A", "A2_A", "HWE_GENO_A", "O(HET)_A", "E(HET)_A", "HWE_P_A", "row")   

hwe_u <- hwe[TEST == 'UNAFF']
hwe_u$row <- as.numeric(rownames(hwe_u))
setkey(hwe_u, row)
names(hwe_u) <- c("CHR_U", "SNP_U", "TEST_U", "A1_U", "A2_U", "HWE_GENO_U", "O(HET)_U", "E(HET)_U", "HWE_P_U", "row")  
hwe_all_a <- hwe_all[hwe_a]
hwe_au <- hwe_all_a[hwe_u]
freq_hwe_au <- freq[hwe_au]
results_freq_hwe_au <- results[freq_hwe_au]
results_freq_hwe_au  %>% filter(SNP != SNP_frq | SNP != SNP_A | SNP != SNP_U )
results_freq_hwe_au  %>%filter(A1 != A1_frq | A1 != A1_A | A1 != A1_U)

rm(freq, hwe, hwe_all, results, hwe_a, hwe_u)


variants <- results_freq_hwe_au[, .(CHR, SNP, BP, A1, A2_frq)]
colnames(variants) <- c("chromosome","snp","GRCh37_bp","A1","A2")
variants$GRCh38_bp <- NA

dbWriteTable(db, name='dimVariant',variants, append = TRUE,row.names = FALSE)

factTable <- results_freq_hwe_au[,.(OR,SE, P, MAF_A, MAF_U, HWE_GENO_A, HWE_GENO_U, HWE_P_A, HWE_P_U, row)]
colnames(factTable) <- c("parameterEstimate","stdError","pValue","MAF_A","MAF_U","HWE_GENO_A", "HWE_GENO_U", "HWE_P_A", "HWE_P_U", 'row')
factTable$parameterType  <-  paramID
factTable$majorAllele  <- 2
factTable$iterationID <- iterationID
factTable$experimentID <- expID
factTable$variantID <- factTable[,row]

dbWriteTable(db, name='factGWAS',factTable[, -'row', with = FALSE], append = TRUE,row.names = FALSE)
