### R code from vignette source 'IntroToAnnotationPackages.Rnw'

###################################################
### code chunk number 1: loadChip
###################################################
library(hgu95av2.db)


###################################################
### code chunk number 2: listContents
###################################################
ls("package:hgu95av2.db")


###################################################
### code chunk number 3: show
###################################################
hgu95av2.db


###################################################
### code chunk number 4: columns
###################################################
columns(hgu95av2.db)


###################################################
### code chunk number 5: help (eval = FALSE)
###################################################
## help("SYMBOL")


###################################################
### code chunk number 6: keytypes
###################################################
keytypes(hgu95av2.db)


###################################################
### code chunk number 7: keys
###################################################
head(keys(hgu95av2.db, keytype="SYMBOL"))


###################################################
### code chunk number 8: selectChip
###################################################
#1st get some example keys
k <- head(keys(hgu95av2.db,keytype="PROBEID"))
# then call select
select(hgu95av2.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="PROBEID")


###################################################
### code chunk number 9: selectOrg1
###################################################
library(org.Hs.eg.db)
columns(org.Hs.eg.db)


###################################################
### code chunk number 10: selectOrg2 (eval = FALSE)
###################################################
## help("SYMBOL") ## for explanation of these columns and keytypes values


###################################################
### code chunk number 11: selectOrg3
###################################################
keytypes(org.Hs.eg.db)
uniKeys <- head(keys(org.Hs.eg.db, keytype="UNIPROT"))
cols <- c("SYMBOL", "PATH")
select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="UNIPROT")


###################################################
### code chunk number 12: selectData
###################################################
load(system.file("extdata", "resultTable.Rda", package="AnnotationDbi"))
head(resultTable)


###################################################
### code chunk number 13: selectOrgData
###################################################
annots <- select(org.Hs.eg.db, keys=rownames(resultTable),
                 columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
resultTable <- merge(resultTable, annots, by.x=0, by.y="ENTREZID")
head(resultTable)


###################################################
### code chunk number 14: selectGO
###################################################
library(GO.db)
GOIDs <- c("GO:0042254","GO:0044183")
select(GO.db, keys=GOIDs, columns="DEFINITION", keytype="GOID")


###################################################
### code chunk number 15: selectTxDb
###################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
columns(txdb)
keytypes(txdb)
keys <- head(keys(txdb, keytype="GENEID"))
cols <- c("TXID", "TXSTART")
select(txdb, keys=keys, columns=cols, keytype="GENEID")



###################################################
### code chunk number 16: SessionInfo
###################################################
sessionInfo()


