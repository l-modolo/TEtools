### R code from vignette source 'create_objects.Rnw'

###################################################
### code chunk number 1: samples
###################################################
tab = data.frame(
  file = c("treated1fb", "treated2fb", "treated3fb", "untreated1fb", "untreated2fb", "untreated3fb", "untreated4fb"),
  type = c("single-read", "paired-end", "paired-end", "single-read", "single-read", "paired-end", "paired-end"),
  "number of lanes" = as.integer(c(5,2,2,2,6,2,2)),
  "total number of reads" = c("35158667", "12242535 (x2)", "12443664 (x2)", "17812866", "34284521", "10542625 (x2)", "12214974 (x 2)"),
  "exon counts" = as.integer(c(15679615, 15620018, 12733865, 14924838, 20764558, 10283129, 11653031)),
  stringsAsFactors = TRUE,
  check.names = FALSE)
library("xtable")
print(xtable(tab, caption = "Read numbers and alignment statistics. The column \\emph{exon counts} refers to the number of reads that could be uniquely aligned to an exon.", label="tab:samples"), file="create_objects_tabsamples.tex")


###################################################
### code chunk number 2: pasilla
###################################################
library("pasilla")


###################################################
### code chunk number 3: extdata
###################################################
inDir = system.file("extdata", package="pasilla", mustWork=TRUE)
dir(inDir)


###################################################
### code chunk number 4: samples1
###################################################
z = sub("fb$", "", tab$file)
nc = nchar(z)
samples = data.frame(
  condition = substr(z, 1, nc-1),
#  replicate = as.integer(substr(z, nc, nc)),
  type = tab$type,
  row.names = tab$file,
  stringsAsFactors = TRUE,
  check.names = FALSE)

#stopifnot(!any(is.na(samples$replicate)))


###################################################
### code chunk number 5: samples
###################################################
head(samples)


###################################################
### code chunk number 6: ecs
###################################################
annotationfile = file.path(inDir, "Dmel.BDGP5.25.62.DEXSeq.chr.gff")


###################################################
### code chunk number 7: read
###################################################
library("DEXSeq")
ecs = read.HTSeqCounts(countfiles = file.path(inDir, paste(rownames(samples), "txt", sep=".")), 
           design = samples, 
           flattenedfile = annotationfile)
sampleNames(ecs) = rownames(samples)


###################################################
### code chunk number 8: gfs
###################################################
genesforsubset = readLines(file.path(inDir, "geneIDsinsubset.txt"))
pasillaExons = subsetByGenes(ecs, genes=genesforsubset)


###################################################
### code chunk number 9: expdata
###################################################
expdata = new("MIAME", 
   name="pasilla knockdown", 
   lab="Genetics and Developmental Biology, University of Connecticut Health Center", 
   contact="Dr. Brenton Graveley", 
   title="modENCODE Drosophila pasilla RNA Binding Protein RNAi knockdown RNA-Seq Studies", 
   url="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE18508", 
   abstract="RNA-seq of 3 biological replicates of from the Drosophila melanogaster 
      S2-DRSC cells that have been RNAi depleted of mRNAs encoding pasilla, a mRNA binding 
      protein and 4 biological replicates of the the untreated cell line.")
   pubMedIds(expdata) <- "20921232"
experimentData(pasillaExons) <- expdata


###################################################
### code chunk number 10: DESeq
###################################################
library("DESeq")
genetable = geneCountTable(ecs)

pasillaGenes = newCountDataSet(genetable, 
  conditions = samples)

experimentData(pasillaGenes) = expdata


###################################################
### code chunk number 11: save (eval = FALSE)
###################################################
## save(pasillaExons, file=file.path("..", "..", "data", "pasillaExons.RData"))
## save(pasillaGenes, file=file.path("..", "..", "data", "pasillaGenes.RData"))


###################################################
### code chunk number 12: sessionInfo
###################################################
toLatex(sessionInfo())


