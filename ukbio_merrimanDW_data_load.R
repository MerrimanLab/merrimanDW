library(RMySQL)
library(data.table)

drv = dbDriver("MySQL")
db = dbConnect(drv, default.file = '~/.my.cnf', dbname="merrimanDW")

dbGetQuery(db,'show tables;')

dbGetQuery(db, 'describe dimSourceDataset;')
dbGetQuery(db, 'select * from dimExperiment')
dbGetQuery(db, 'select * from dimIteration')

i = 1
iterationID = i
datasetID <- dbGetQuery(db, 'select datasetID FROM dimSourceDataset WHERE datasetName = "UKBIOBANK";')

param = "OR"
paramID <-  dbGetQuery(db,paste0('select parameterTypeID FROM dimParameter WHERE parameterTypeDesc = ', '"', param,'"', ';'))

caseCond <- "All"
controlCond <- "Not diuretic"
expID <- dbGetQuery(db, paste0('select experimentID from dimExperiment WHERE caseCondition = ', '"',caseCond, '"','and controlCondition = ', '"',controlCond, '"' ,' and datasetID =', datasetID,';'))


dbGetQuery(db, paste0('INSERT into dimIteration (experimentID, iterationID) values (',expID ,',' , iterationID,');'))

fam_file <- fread('/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/controls_nd/all/chrallimpv1.fam_goutall_nd1', header=FALSE)
colnames(fam_file) <- c('fid','iid','PID','MID','SEX','AFF')

fam_file[AFF == 1, .(fid, iid)]
fam_file[AFF == 2, .(fid,iid)]
fam_file$experimentID <-  expID
fam_file$iterationID <- iterationID

dbWriteTable(db, name = 'dimPeopleCase', fam_file[AFF == 2, .(fid,iid, experimentID)], append=TRUE,row.names = FALSE)
dbWriteTable(db, name = 'dimPeopleControl', fam_file[AFF == 1, .(fid,iid, experimentID, iterationID)], append = TRUE ,row.names = FALSE)

results <- fread('/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/controls_nd/all/results/goutall1_chr1.assoc.logistic', header=TRUE)
colnames(results)
results <- results[TEST == 'ADD']
results$row <- as.numeric(rownames(results))
setkey(results, row)


freq <- fread('/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/controls_nd/all/results/goutall1_chr1.frq.cc', header=TRUE)
freq$row <- as.numeric(rownames(freq))
setkey(freq, row)
colnames(freq) <- c("CHR_frq","SNP_frq","A1_frq","A2_frq","MAF_A","MAF_U","NCHROBS_A","NCHROBS_U","row")      

hwe <- fread('/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/controls_nd/all/results/goutall1_chr1.hwe', header=TRUE)
hwe_all <- hwe[TEST == 'ALL']
hwe_all$row <- as.numeric(rownames(hwe_all))
setkey(hwe_all, row)


hwe_a <- hwe[TEST == 'AFF']
hwe_a$row <- as.numeric(rownames(hwe_a))
setkey(hwe_a, row)
names(hwe_a) <- c("CHR_A", "SNP_A", "TEST_A", "A1_A", "A2_A", "HWE_GENO_A", "O(HET)_A", "E(HET)_A", "HWE_P_A", "row")   

hwe_u <- hwe[TEST == 'UNAFF']
hwe_u$row <- as.numeric(rownames(hwe_u))
setkey(hwe_u, row)
names(hwe_u) <- c("CHR_U", "SNP_U", "TEST_U", "A1_U", "A2_U", "HWE_GENO_U", "O(HET)_U", "E(HET)_U", "HWE_P_U", "row")  

hwe_au <- hwe_a[hwe_u]
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
