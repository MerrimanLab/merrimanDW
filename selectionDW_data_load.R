library(RMySQL)

drv = dbDriver("MySQL")
selection_db = dbConnect(drv,default.file = '~/.my.cnf', dbname="selectionDW")
scratch_dir <- '/media/xsan/scratch/merrimanlab/murray/working_dir/'


###
### Populations
###
populations <- read.csv('populations.csv', header=TRUE)
populations$popID <- as.numeric(row.names(populations))
dbWriteTable(selection_db, name = 'dimPopData', populations, row.names = FALSE, append = TRUE)

pops <- dbGetQuery(selection_db, "select * from dimPopData")
pop_1kg <- pops[1:26,]

###
### Experiments
###
dbGetQuery(selection_db, 'insert into dimExperiment (description) VALUES ("AXIOM_unimputed"),("OMNI_unimputed"), ("AXIOM and OMNI imputed info score threshold 0.3"), ("AXIOM and OMNI imputed info score threshold 0.8"), ("coreExome_v24_unimputed")')

###
### Genes
###
genes <- read.csv('dimGene.csv', header=TRUE)
var <- genes[, c('chromosome_name', 'start_position','end_position')]
names(var) <- c('chrom','chrom_start','chrom_end')
var$posID <- as.numeric(row.names(var))
dbWriteTable(selection_db, name = 'dimPos', var, append=TRUE, row.names = FALSE)
gen <- genes[,c('ensembl_gene_id','hgnc_symbol')]
names(gen) <- c('EnsGeneID', "GeneName")
gen$posID <- as.numeric(row.names(gen))
dbWriteTable(selection_db, name = 'dimGene', gen, append = TRUE, row.names = FALSE)


###
### STATS
###
# Tajimas D
dbGetQuery(selection_db, 'insert into dimStat (statName, statDescription) values ("TajimasD", "Tajimas D"), ("NumSites_TajimasD", "Number of sites in window for Tajimas D")')
dbGetQuery(selection_db, 'insert into dimStat (statName, statDescription) values ("NumSites_FayWuH","Number of sites in window for Fay and Wus H"),("s","Segregating sites"), ("Eta","NULL"),("Eta_e","NULL"),("Pi","NULL"),("FuLiD","Fu and Lis D"),("FuLiF","Fu and Lis F"),("FayWuH","Fay and Wus H") ')


#### tajimas d - axiom unimputed ####
##
exID <- 1
dimStat <- dbReadTable(selection_db, name = 'dimStat')
td <- which(dimStat$statName == "TajimasD")
nsite <- which(dimStat$statName == 'NumSites_TajimasD')
tajd_path <- paste0(scratch_dir,"/nesi_retrieved/Unimputed/Axiom/TD/")
for (pop in pops$pop[1:27]){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(tajd_path,'/',pop,'/results/',pop,chr,".taj_d"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
    tmp$popID <- as.numeric(which(pops$pop==pop))
    tmp
    #dimVar
    tmp$chrom <- chr
    tmp$chrom_start <- tmp$BIN_START
    w <- tmp$BIN_START[3] - tmp$BIN_START[2]
    tmp$chrom_end <- tmp$BIN_START +w -1
    tmp$slide = w
    var <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(var, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    newVar$slide <- newVar$chrom_end - newVar$chrom_start +1
    newVar <- newVar[is.na(newVar$posID),]
    if(nrow(newVar) > 0){
    dbWriteTable(selection_db, name = 'dimPos', newVar[, c('chrom_start','chrom_end','chrom', 'slide')], append = TRUE, row.names=FALSE)
    }
    #update posID
    updatedVar <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(updatedVar, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    # insert tajimas d
    newVar$experimentID <- exID
    newVar$statValue <- newVar$TajimaD
    newVar$statID <- td
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    # insert number of sites
    newVar$statValue <- newVar$N_SNPS
    newVar$statID <- nsite
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
                                                    
  }
}

#### tajimas d - omni unimputed ####

exID <- 2
dimStat <- dbReadTable(selection_db, name = 'dimStat')
td <- which(dimStat$statName == "TajimasD")
nsite <- which(dimStat$statName == 'NumSites_TajimasD')
tajd_path <- paste0(scratch_dir,"/nesi_retrieved/Unimputed/Omni/TD/")
for (pop in pops$pop[c(1:26,28)]){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(tajd_path,'/',pop,'/results/',pop,chr,".taj_d"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
    tmp$popID <- as.numeric(which(pops$pop==pop))
    tmp
    #dimVar
    tmp$chrom <- chr
    tmp$chrom_start <- tmp$BIN_START
    w <- tmp$BIN_START[3] - tmp$BIN_START[2]
    tmp$chrom_end <- tmp$BIN_START +w -1
    tmp$slide = w
    var <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(var, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    newVar$slide <- newVar$chrom_end - newVar$chrom_start +1
    newVar <- newVar[is.na(newVar$posID),]
    if(nrow(newVar) > 0){
      dbWriteTable(selection_db, name = 'dimPos', newVar[, c('chrom_start','chrom_end','chrom', 'slide')], append = TRUE, row.names=FALSE)
    }
    #update posID
    updatedVar <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(updatedVar, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    # insert tajimas d
    newVar$experimentID <- exID
    newVar$statValue <- newVar$TajimaD
    newVar$statID <- td
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    # insert number of sites
    newVar$statValue <- newVar$N_SNPS
    newVar$statID <- nsite
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
  }
}


#### tajimas d - 1KG pops - axiom_omni info 0.3 ####

exID <- NA
dimStat <- dbReadTable(selection_db, name = 'dimStat')
td <- which(dimStat$statName == "TajimasD")
nsite <- which(dimStat$statName == 'NumSites_TajimasD')

tajd_path <- "/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/Phase3_selection_results/TajimaD/"
for (pop in pops$pop[1:26]){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(tajd_path,'/',pop,'/results/',pop,chr,".taj_d"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
    tmp$popID <- as.numeric(which(pops$pop==pop))
    tmp
    #dimVar
    tmp$chrom <- chr
    tmp$chrom_start <- tmp$BIN_START
    w <- tmp$BIN_START[3] - tmp$BIN_START[2]
    tmp$chrom_end <- tmp$BIN_START +w -1
    tmp$slide = w
    var <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(var, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    newVar$slide <- newVar$chrom_end - newVar$chrom_start +1
    newVar <- newVar[is.na(newVar$posID),]
    if(nrow(newVar) > 0){
      dbWriteTable(selection_db, name = 'dimPos', newVar[, c('chrom_start','chrom_end','chrom', 'slide')], append = TRUE, row.names=FALSE)
    }
    #update posID
    updatedVar <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(updatedVar, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    # insert tajimas d
    newVar$experimentID <- exID
    newVar$statValue <- newVar$TajimaD
    newVar$statID <- td
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    # insert number of sites
    newVar$statValue <- newVar$N_SNPS
    newVar$statID <- nsite
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
  }
}


#### tajimas d - axiom_omni info 0.3 ####

exID <- NA
dimStat <- dbReadTable(selection_db, name = 'dimStat')
td <- which(dimStat$statName == "TajimasD")
nsite <- which(dimStat$statName == 'NumSites_TajimasD')
tajd_path <- paste0(scratch_dir,"/nesi_retrieved/Unimputed/Omni/TD/")
for (pop in pops$pop[27]){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(tajd_path,'/',pop,'/results/',pop,chr,".taj_d"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
    tmp$popID <- as.numeric(which(pops$pop==pop))
    tmp
    #dimVar
    tmp$chrom <- chr
    tmp$chrom_start <- tmp$BIN_START
    w <- tmp$BIN_START[3] - tmp$BIN_START[2]
    tmp$chrom_end <- tmp$BIN_START +w -1
    tmp$slide = w
    var <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(var, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    newVar$slide <- newVar$chrom_end - newVar$chrom_start +1
    newVar <- newVar[is.na(newVar$posID),]
    if(nrow(newVar) > 0){
      dbWriteTable(selection_db, name = 'dimPos', newVar[, c('chrom_start','chrom_end','chrom', 'slide')], append = TRUE, row.names=FALSE)
    }
    #update posID
    updatedVar <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(updatedVar, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    # insert tajimas d
    newVar$experimentID <- exID
    newVar$statValue <- newVar$TajimaD
    newVar$statID <- td
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    # insert number of sites
    newVar$statValue <- newVar$N_SNPS
    newVar$statID <- nsite
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
  }
}


#### Fay Wu's H - Axiom unimputed ####
exID <- 1
faw_path <- paste0(scratch_dir,"/nesi_retrieved/Unimputed/Axiom/FAWH/")
dimStat <- dbReadTable(selection_db, name = 'dimStat')
s <- which(dimStat$statName == "s")
nsite <- which(dimStat$statName == 'NumSites_FayWuH')
eta <- which(dimStat$statName == "Eta")
eta_e <- which(dimStat$statName == "Eta_e")
pi <- which(dimStat$statName == "Pi")
fuli_d <-which(dimStat$statName == "FuLiD")
fuli_f <-which(dimStat$statName == "FuLiF")
faywuh <- which(dimStat$statName == "FayWuH")
for (pop in pops$pop[c(1:27)]){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(faw_path,'/',pop,'/results/',pop,chr,".faw"), header=FALSE, stringsAsFactors = FALSE, skip =5)
    names(tmp)=c("RefStart","Refend","RefMid","chrom_start","chrom_end","Midpoint","num_sites","missing","s","eta","eta_e","pi","fuli_d","fuli_f","faywu_h")
    tmp$popID <- which(pops$pop==pop)
    tmp$chrom <- chr
    w <- tmp$chrom_end[3] - tmp$chrom_start[3] +1
    tmp$slide = w
    
    var <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(var, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    newVar$slide <- newVar$chrom_end - newVar$chrom_start +1
    dbWriteTable(selection_db, name = 'dimPos', newVar[, c('chrom_start','chrom_end','chrom', 'slide')], append = TRUE, row.names=FALSE)
    #update posID
    updatedVar <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(updatedVar, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)

    # insert S
    newVar$experimentID <- exID
    newVar$statValue <- newVar$s
    newVar$statID <- s
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    # insert number of sites
    newVar$statValue <- newVar$num_sites
    newVar$statID <- nsite
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)

    # insert eta
    newVar$statValue <- newVar$num_sites
    newVar$statID <- eta
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert eta_e
    newVar$statValue <- newVar$eta_e
    newVar$statID <- eta_e
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert pi
    newVar$statValue <- newVar$pi
    newVar$statID <- pi
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert fu li's d
    newVar$statValue <- newVar$fuli_d
    newVar$statID <- fuli_d
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert fu li's f
    newVar$statValue <- newVar$fuli_f
    newVar$statID <- fuli_f
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert fay wu's h
    newVar$statValue <- newVar$faywu_h
    newVar$statID <- faywuh
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
  }
}


#### Fay Wu's H - OMNI unimputed ####
exID <- 1
faw_path <- paste0(scratch_dir,"/nesi_retrieved/Unimputed/Omni/FAWH/")
dimStat <- dbReadTable(selection_db, name = 'dimStat')
s <- which(dimStat$statName == "s")
nsite <- which(dimStat$statName == 'NumSites_FayWuH')
eta <- which(dimStat$statName == "Eta")
eta_e <- which(dimStat$statName == "Eta_e")
pi <- which(dimStat$statName == "Pi")
fuli_d <-which(dimStat$statName == "FuLiD")
fuli_f <-which(dimStat$statName == "FuLiF")
faywuh <- which(dimStat$statName == "FayWuH")
for (pop in pops$pop[c(1:27)]){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(faw_path,'/',pop,'/results/',pop,chr,".faw"), header=FALSE, stringsAsFactors = FALSE, skip =5)
    names(tmp)=c("RefStart","Refend","RefMid","chrom_start","chrom_end","Midpoint","num_sites","missing","s","eta","eta_e","pi","fuli_d","fuli_f","faywu_h")
    tmp$popID <- which(pops$pop==pop)
    tmp$chrom <- chr
    w <- tmp$chrom_end[3] - tmp$chrom_start[3] +1
    tmp$slide = w
    
    var <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(var, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    newVar$slide <- newVar$chrom_end - newVar$chrom_start
    dbWriteTable(selection_db, name = 'dimPos', newVar[, c('chrom_start','chrom_end','chrom', 'slide')], append = TRUE, row.names=FALSE)
    #update posID
    updatedVar <- dbReadTable(selection_db, name = 'dimPos')
    newVar <- merge(updatedVar, tmp, by = c('chrom','chrom_start','chrom_end', 'slide'), all.y = TRUE)
    
    # insert S
    newVar$experimentID <- exID
    newVar$statValue <- newVar$s
    newVar$statID <- s
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    # insert number of sites
    newVar$statValue <- newVar$num_sites
    newVar$statID <- nsite
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert eta
    newVar$statValue <- newVar$num_sites
    newVar$statID <- eta
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert eta_e
    newVar$statValue <- newVar$eta_e
    newVar$statID <- eta_e
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert pi
    newVar$statValue <- newVar$pi
    newVar$statID <- pi
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert fu li's d
    newVar$statValue <- newVar$fuli_d
    newVar$statID <- fuli_d
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert fu li's f
    newVar$statValue <- newVar$fuli_f
    newVar$statID <- fuli_f
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
    
    # insert fay wu's h
    newVar$statValue <- newVar$faywu_h
    newVar$statID <- faywuh
    dbWriteTable(selection_db, name = 'intraSel', newVar[,c('popID', 'posID','statValue','statID','experimentID' )], row.names=FALSE, append = TRUE)
  }
}


# daf
# DAF
# <pos> <Ref> <Alt> <Anc> <MAF> <DAF>


exID <- 1
dimStat <- dbReadTable(selection_db, name = 'dimStat')

daf_path <- paste0(scratch_dir,"/nesi_retrieved/Unimputed/Axiom/AF/")
for (pop in pops$code){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(daf_path,"/",pop,'/',pop,chr,"_aachanged.af"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
    tmp$pop <- which(pops$code==pop)
    tmp$chrom <- chr
    tmp$chrom_start <- tmp$Pos
    tmp$chrom_end <- tmp$Pos +1
    
    dbWriteTable(axiom_db, name = "allele_freq", tmp[,c("chrom","chrom_start","chrom_end", "Ref", "Alt", "Anc", "MAF","DAF","pop")], overwrite = FALSE, append=TRUE, row.names=FALSE)
  }
}



##
# KaKs
##
kaks_path <- "/media/scratch/merrimanlab/murray/working_dir/Phase3_selection_results/KaKs/" 
for (pop in pop_1kg[,2]){
  for(chr in 1:22){
    tmp <- read.table(paste0(kaks_path,pop,chr,".kaks"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
    names(tmp) <- c("gene_id", "gene_name", "ka", "ks", "ka_div_ks_1")
    tmp$pop <- which(pops$code==pop)
    tmp$chrom <- chr
    
    dbWriteTable(db, name = "kaks", tmp[,c("chrom","gene_id","gene_name","ka","ks" ,"pop")], overwrite = FALSE, append=TRUE, row.names=FALSE)
  }
}
for (pop in c("AXIOM","OMNI")){
  for(info in c(0.3,0.8)){
    for(chr in 1:22){
      tmp <- read.table(paste0(kaks_path,"Info",info,"/",pop,chr,".kaks"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
      names(tmp) <- c("gene_id", "gene_name", "ka", "ks", "ka_div_ks_1")
      tmp$pop <- which(pops$code==paste0(pop,info))
      tmp$chrom <- chr
      
      
      dbWriteTable(db, name = "kaks", tmp[,c("chrom","gene_id","gene_name","ka","ks" ,"pop")], overwrite = FALSE, append=TRUE, row.names=FALSE)
    }
  }
}



# dbGetQuery(db, "CREATE TABLE `faw`(`chrom` smallint(3) unsigned,
#            `chrom_start` int(10) unsigned,
#            `chrom_end` int(10) unsigned,
#             `num_sites` int (10),
#             `missing` int(10),
#            `s` int(10),
#            `eta` int(10),
#            `eta_e` int(10), 
#            `pi` float,
#            `fuli_d` float,
#            `fuli_f` float,
#            `faywu_h` float,
#            `pop` smallint(3) unsigned,
#            `window_id` smallint(3) unsigned,
#             FOREIGN KEY (pop)
#            references population (id),
#           
#             FOREIGN KEY (window_id)
#               references window_info (id)
#            ) ENGINE = INNODB;")
faw_path <- "/media/scratch/merrimanlab/murray/working_dir/nesi_retrieved/Unimputed/Axiom/FAWH/"
for (pop in pops$code){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(faw_path,'/',pop,'/results/',pop,chr,".faw"), header=FALSE, stringsAsFactors = FALSE, skip =5)
    names(tmp)=c("RefStart","Refend","RefMid","chrom_start","chrom_end","Midpoint","num_sites","missing","s","eta","eta_e","pi","fuli_d","fuli_f","faywu_h")
    tmp$pop <- which(pops$code==pop)
    tmp$chrom <- chr
    
    w <- tmp$chrom_end[3] - tmp$chrom_start[3] +1
    tmp$window_id <- window_info[which(window_info$width == w & window_info$slide == w),'id']
    tmp <- replace(tmp, is.na(tmp), NA)
    dbWriteTable(axiom_db, name = "faw", tmp[,c("chrom","chrom_start","chrom_end", "num_sites","missing", "s", "eta", "eta_e","pi","fuli_d","fuli_f","faywu_h" ,"pop", "window_id")], overwrite = FALSE, append=TRUE, row.names=FALSE)
  }
}


#ihs
#<locusID> <physicalPos> <'1' freq> <ihh1> <ihh0> <unstandardized iHS> <norm> <crit>
ihs_path <- "/media/scratch/merrimanlab/murray/working_dir/nesi_retrieved/Unimputed/Axiom/iHS/"
for (pop in pops$code){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(ihs_path,pop,"_axiom_",chr,".ihs.out.100bins.norm"), header=FALSE, stringsAsFactors = FALSE, sep="\t")
    names(tmp) <- c("locus_id", "chrom_start", "freq_1", "ihh1","ihh0", "unstd_ihs", "norm_ihs", "significant")
    tmp$pop <- which(pops$code==pop)
    tmp$chrom <- chr
    
    tmp$chrom_end <- tmp$chrom_start +1
    dbWriteTable(axiom_db, name = "ihs", tmp[,c("chrom","chrom_start","chrom_end", "locus_id", "freq_1", "ihh0", "ihh1", "unstd_ihs", "norm_ihs", "significant" ,"pop")], overwrite = FALSE, append=TRUE, row.names=FALSE)
  }
}




# nSL
# dbGetQuery(db, "CREATE TABLE `nsl`(`chrom` smallint(3) unsigned,
#            `chrom_start` int(10) unsigned,
#            `chrom_end` int(10) unsigned,
#             `locus_id` varchar (20),
#            `freq_1` float,
#             `sL1` float,
#             `sL0` float,
#             `unstd_nsl` float,
#             `norm_nsl` float,
#             `significant` boolean, 
#            `pop` smallint(3) unsigned, 
# 
#             FOREIGN KEY (pop)
#               references population (id)
#             )ENGINE = INNODB;")
# <locusID> <physicalPos> <'1' freq> <ihh1> <ihh0> <unstandardized iHS> <norm> <crit>
nsl_path <- "/media/scratch/merrimanlab/murray/working_dir/nesi_retrieved/Unimputed/Axiom/nSL/"
for (pop in pops$code){
  for(chr in 1:22){
    print(paste0(pop,chr))
    tmp <- read.table(paste0(nsl_path,pop,"_axiom_",chr,".nsl.out.100bins.norm"), header=FALSE, stringsAsFactors = FALSE, sep="\t")
    names(tmp) <- c("locus_id", "chrom_start", "freq_1", "sL1","sL0", "unstd_nsl", "norm_nsl", "significant")
    tmp$pop <- which(pops$code==pop)
    tmp$chrom <- chr
    
    tmp$chrom_end <- tmp$chrom_start +1
    dbWriteTable(axiom_db, name = "nsl", tmp[,c("chrom","chrom_start","chrom_end", "locus_id", "freq_1", "sL1", "sL0", "unstd_nsl", "norm_nsl", "significant" ,"pop")], overwrite = FALSE, append=TRUE, row.names=FALSE)
  }
}





# dbGetQuery(db, "CREATE TABLE `xpehh_axiom`(`chrom` smallint(3) unsigned,
#            `chrom_start` int(10) unsigned,
#            `chrom_end` int(10) unsigned,
#             `gpos` float,
#            `popA_1_freq` float,
#            `ihhA` float,
#            `popB_1_freq` float,
#            `ihhB` float,
#            `unstd_xpehh` float,
#            `norm_xpehh` float,
#            `significant` boolean, 
#            `popA` smallint(3) unsigned, 
#            `popB` smallint(3) unsigned,
#            
#            FOREIGN KEY (popA)
#            references population (id),
#            FOREIGN KEY (popB)
#            references population (id)
# )ENGINE = INNODB;")
# <locusID> <physicalPos> <geneticPos> <popA '1' freq> <ihhA> <popB '1' freq> <ihhB> <unstandardized XPEHH> <norm XPEHH> <significant (0/1)>
# chr  pos     gpos    p1      ihh1    p2      ihh2    xpehh   normxpehh       crit
xpehh_path <- ""
for (pop2 in pop_1kg[,2]){
  pop = "AXIOM"
  for(info in c(0.3,0.8)){
    for(chr in 1:22){
      tmp <- read.table(paste0(xpehh_path,"Info",info,"/",pop,"_",pop2,'_',chr,".xpehh.out.norm"), header=FALSE, stringsAsFactors = FALSE, sep="\t", skip=1)
      names(tmp) <- c('chrom', "chrom_start",'gpos', "popA_1_freq", "ihhA",'popB_1_freq',"ihhB", "unstd_xpehh", "norm_xpehh", "significant")
      tmp$popA <- which(pops$code==paste0(pop,info))
      tmp$popB <- which(pops$code == pop2)
      
      tmp$chrom_end <- tmp$chrom_start +1
      dbWriteTable(db, name = "xpehh_axiom", tmp[,c("chrom","chrom_start","chrom_end",'gpos', "popA_1_freq", "popB_1_freq","ihhA", "ihhB", "unstd_xpehh", "norm_xpehh", "significant" ,"popA","popB")], overwrite = FALSE, append=TRUE, row.names=FALSE)
    }
  }
}

pop="AXIOM"
pop2="OMNI"
for(info in c(0.3,0.8)){
  for(chr in 1:22){
    tmp <- read.table(paste0(xpehh_path,"Info",info,"/",pop,"_",pop2,'_',chr,".xpehh.out.norm"), header=FALSE, stringsAsFactors = FALSE, sep="\t", skip=1)
    names(tmp) <- c('chrom', "chrom_start",'gpos', "popA_1_freq", "ihhA",'popB_1_freq',"ihhB", "unstd_xpehh", "norm_xpehh", "significant")
    tmp$popA <- which(pops$code==paste0(pop,info))
    tmp$popB <- which(pops$code == paste0(pop2,info))
    
    tmp$chrom_end <- tmp$chrom_start +1
    dbWriteTable(db, name = "xpehh_axiom", tmp[,c("chrom","chrom_start","chrom_end",'gpos', "popA_1_freq", "popB_1_freq","ihhA", "ihhB", "unstd_xpehh", "norm_xpehh", "significant" ,"popA","popB")], overwrite = FALSE, append=TRUE, row.names=FALSE)
  }
}



# dbGetQuery(db, "CREATE TABLE `xpehh_omni`(`chrom` smallint(3) unsigned,
#            `chrom_start` int(10) unsigned,
#            `chrom_end` int(10) unsigned,
#            `gpos` float,
#            `popA_1_freq` float,
#            `ihhA` float,
#            `popB_1_freq` float,
#            `ihhB` float,
#            `unstd_xpehh` float,
#            `norm_xpehh` float,
#            `significant` boolean, 
#            `popA` smallint(3) unsigned, 
#            `popB` smallint(3) unsigned,
#            
#            FOREIGN KEY (popA)
#            references population (id),
#            FOREIGN KEY (popB)
#            references population (id)
# )ENGINE = INNODB;")
#<locusID> <physicalPos> <geneticPos> <popA '1' freq> <ihhA> <popB '1' freq> <ihhB> <unstandardized XPEHH> <norm XPEHH> <significant (0/1)>
# chr  pos     gpos    p1      ihh1    p2      ihh2    xpehh   normxpehh       crit
xpehh_path <- ""
for (pop2 in pop_1kg[,2]){
  pop = "OMNI"
  for(info in c(0.3,0.8)){
    for(chr in 1:22){
      tmp <- read.table(paste0(xpehh_path,"Info",info,"/",pop,"_",pop2,'_',chr,".xpehh.out.norm"), header=TRUE, stringsAsFactors = FALSE, sep="\t")
      names(tmp) <- c('chrom', "chrom_start",'gpos', "popA_1_freq", "ihhA",'popB_1_freq',"ihhB", "unstd_xpehh", "norm_xpehh", "significant")
      tmp$popA <- which(pops$code==paste0(pop,info))
      tmp$popB <- which(pops$code == pop2)
      
      tmp$chrom_end <- tmp$chrom_start +1
      dbWriteTable(db, name = "xpehh_omni", tmp[,c("chrom","chrom_start","chrom_end",'gpos', "popA_1_freq", "popB_1_freq","ihhA", "ihhB", "unstd_xpehh", "norm_xpehh", "significant" ,"popA","popB")], overwrite = FALSE, append=TRUE, row.names=FALSE)
    }
  }
}




##
## Fst
##
CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST
dbGetQuery(db, "CREATE TABLE `fst_axiom`(`chrom` smallint(3) unsigned,
           `chrom_start` int(10) unsigned,
           `chrom_end` int(10) unsigned,
           `num_snps` int (10),
           `weighted_fst` float,
           `mean_fst` float,
           `window_id` smallint(3) unsigned,
           `popA` smallint(3) unsigned, 
           `popB` smallint(3) unsigned, 
           
           FOREIGN KEY (popA)
           references population (id),
           FOREIGN KEY (popB)
           references population (id),
           FOREIGN KEY (window_id)
           references window_info (id)
)ENGINE = INNODB;")

fst_path <- "~/Murray/Bioinformatics/working_dir/Phase3_selection_results/FST/"
# CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST
window_info <- dbGetQuery(db, 'select * from window_info')
for (pop2 in pop_1kg[,2]){
  pop = "AXIOM"
  for(info in c(0.3,0.8)){
    print(pop2)
    for(chr in 1:22){
      print(chr)
      tmp <- read.table(paste0(fst_path,"Info",info,"/",pop,"_",pop2,chr,".windowed.weir.fst"), header=FALSE, stringsAsFactors = FALSE, sep="\t", skip=1)
      names(tmp) <- c("chrom", "chrom_start", "chrom_end", "num_snps", "weighted_fst", "mean_fst")
      w = tmp$chrom_end[3] - tmp$chrom_start[3] +1
      s = tmp$chrom_start[4] - tmp$chrom_start[3]
      tmp$window_id <- window_info[which(window_info$width == w & window_info$slide == s),'id']
      tmp$popA <- which(pops$code==paste0(pop,info))
      tmp$popB <- which(pops$code == pop2)
      
      dbWriteTable(db, name = "fst_axiom", tmp[,c("chrom","chrom_start","chrom_end","num_snps", "weighted_fst", "mean_fst", "window_id" ,"popA","popB")], overwrite = FALSE, append=TRUE, row.names=FALSE)
    }
  }
}

for (pop2 in pop_1kg[,2]){
  pop = "OMNI"
  for(info in c(0.3,0.8)){
    for(chr in 1:22){
      tmp <- read.table(paste0(fst_path,"Info",info,"/",pop,"_",pop2,chr,".windowed.weir.fst"), header=FALSE, stringsAsFactors = FALSE, sep="\t", skip=1)
      names(tmp) <- c("chrom", "chrom_start", "chrom_end", "num_snps", "weighted_fst", "mean_fst")
      w = tmp$chrom_end[3] - tmp$chrom_start[3] +1
      s = tmp$chrom_start[4] - tmp$chrom_start[3]
      tmp$window_id <- window_info[which(window_info$width == w & window_info$slide == s),'id']
      tmp$popA <- which(pops$code==paste0(pop,info))
      tmp$popB <- which(pops$code == pop2)
      
      dbWriteTable(db, name = "fst_omni", tmp[,c("chrom","chrom_start","chrom_end","num_snps", "weighted_fst", "mean_fst", "window_id" ,"popA","popB")], overwrite = FALSE, append=TRUE, row.names=FALSE)
    }
  }
}



for (pop2 in c("OMNI")){
  pop = "AXIOM"
  for(info in c(0.3,0.8)){
    for(chr in 1:22){
      tmp <- read.table(paste0(fst_path,"Info",info,"/",pop,"_",pop2,chr,".windowed.weir.fst"), header=FALSE, stringsAsFactors = FALSE, sep="\t", skip=1)
      names(tmp) <- c("chrom", "chrom_start", "chrom_end", "num_snps", "weighted_fst", "mean_fst")
      w = tmp$chrom_end[3] - tmp$chrom_start[3] +1
      s = tmp$chrom_start[4] - tmp$chrom_start[3]
      tmp$window_id <- window_info[which(window_info$width == w & window_info$slide == s),'id']
      tmp$popA <- which(pops$code==paste0(pop,info))
      tmp$popB <- which(pops$code == paste0(pop2,info))
      
      dbWriteTable(db, name = "fst_axiom", tmp[,c("chrom","chrom_start","chrom_end","num_snps", "weighted_fst", "mean_fst", "window_id" ,"popA","popB")], overwrite = FALSE, append=TRUE, row.names=FALSE)
    }
  }
}


for (pop2 in c("AXIOM")){
  pop = "OMNI"
  for(info in c(0.3,0.8)){
    for(chr in 1:22){
      tmp <- read.table(paste0(fst_path,"Info",info,"/",pop,"_",pop2,chr,".windowed.weir.fst"), header=FALSE, stringsAsFactors = FALSE, sep="\t", skip=1)
      names(tmp) <- c("chrom", "chrom_start", "chrom_end", "num_snps", "weighted_fst", "mean_fst")
      w = tmp$chrom_end[3] - tmp$chrom_start[3] +1
      s = tmp$chrom_start[4] - tmp$chrom_start[3]
      tmp$window_id <- window_info[which(window_info$width == w & window_info$slide == s),'id']
      tmp$popA <- which(pops$code==paste0(pop,info))
      tmp$popB <- which(pops$code == paste0(pop2,info))
      
      dbWriteTable(db, name = "fst_omni", tmp[,c("chrom","chrom_start","chrom_end","num_snps", "weighted_fst", "mean_fst", "window_id" ,"popA","popB")], overwrite = FALSE, append=TRUE, row.names=FALSE)
    }
  }
}