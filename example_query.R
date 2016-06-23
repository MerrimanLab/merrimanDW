library(RMySQL)
library(data.table)
library(dplyr)
drv = dbDriver("MySQL")
db = dbConnect(drv, default.file = '~/.my.cnf', dbname="merrimanDW_test")

dbGetQuery()

ukvariant <- dbGetQuery(db, paste('select * from dimVariant where snp in (', paste(paste0("'",kot_t1$SNP,"'"), collapse = ','),")"))



dbGetQuery(db, "select * from dimExperiment where caseCondition IN (1,6) AND controlCondition = 1")
re <- dbGetQuery(db, "select  dV.snp, dV.chromosome, dV.GRCh37_bp, dV.A1, dV.A2, fG.*, dE.Description from factGWAS as fG join dimVariant as dV on fG.variantID = dV.variantID join dimExperiment as dE on dE.experimentID = fG.experimentID where dV.chromosome = 3 AND dV.GRCh37_bp BETWEEN 195776000 and 195835000 AND dE.caseCondition IN (1,6) AND dE.controlCondition = 1")

re <- re[, !(names(re) %in% c("variantID", "infoscore","experimentID"))]
write.table(re, file = '~/TFRC.csv', col.names = TRUE, row.names = FALSE, quote = FALSE, sep=",")
