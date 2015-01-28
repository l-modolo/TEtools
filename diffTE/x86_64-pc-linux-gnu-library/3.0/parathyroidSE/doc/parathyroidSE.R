### R code from vignette source 'parathyroidSE.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: options
###################################################
options(digits=3, width=80, prompt=" ", continue=" ")


###################################################
### code chunk number 3: getExonsByGene (eval = FALSE)
###################################################
## library("GenomicFeatures")
## hse <- makeTranscriptDbFromBiomart(biomart="ensembl",
##                                    dataset="hsapiens_gene_ensembl")
## exonsByGene <- exonsBy(hse, by="gene")


###################################################
### code chunk number 4: loadExonsByGene
###################################################
library("parathyroidSE")
data(exonsByGene)


###################################################
### code chunk number 5: locateFiles
###################################################
bamDir <- system.file("extdata",package="parathyroidSE",mustWork=TRUE)
fls <- list.files(bamDir, pattern="bam$",full=TRUE)


###################################################
### code chunk number 6: bamFiles
###################################################
library("Rsamtools")
bamLst <- BamFileList(fls, index=character(), 
                      yieldSize=100000, obeyQname=TRUE)


###################################################
### code chunk number 7: countGenes
###################################################
parathyroidGenesSE <- summarizeOverlaps(exonsByGene, bamLst, 
                                        mode="Union", 
                                        singleEnd=FALSE, 
                                        ignore.strand=TRUE,
                                        fragments=TRUE)


###################################################
### code chunk number 8: getExonsByTranscript (eval = FALSE)
###################################################
## exonicParts <- disjointExons(hse)


###################################################
### code chunk number 9: importExonicParts
###################################################
data(exonicParts)


###################################################
### code chunk number 10: exonicPartsResult
###################################################
exonicParts[1:3]


###################################################
### code chunk number 11: exonCounts
###################################################
parathyroidExonsSE <- summarizeOverlaps(exonicParts, bamLst,
                                        mode="Union",
                                        singleEnd=FALSE,
                                        ignore.strand=TRUE,
                                        inter.feature=FALSE,
                                        fragments=TRUE)


###################################################
### code chunk number 12: metaExons
###################################################
str(metadata(rowData(parathyroidGenesSE)))


###################################################
### code chunk number 13: getGEO
###################################################
library("GEOquery")
gse37211 <- getGEO(filename=system.file("extdata/GSE37211_series_matrix.txt",
                               package="parathyroidSE",mustWork=TRUE))
samples <- pData(gse37211)[,c("characteristics_ch1","characteristics_ch1.2",
                              "characteristics_ch1.3","relation")]
colnames(samples) <- c("patient","treatment","time","experiment")
samples$patient <- sub("patient: (.+)","\\1",samples$patient)
samples$treatment <- sub("agent: (.+)","\\1",samples$treatment)
samples$time <- sub("time: (.+)","\\1",samples$time)
samples$experiment <- sub("SRA: http://www.ncbi.nlm.nih.gov/sra\\?term=(.+)","\\1",
                          samples$experiment)
samples


###################################################
### code chunk number 14: getSRA (eval = FALSE)
###################################################
## library("SRAdb")
## sqlfile <- getSRAdbFile()
## sra_con <- dbConnect(SQLite(),sqlfile)
## conversion <- sraConvert(in_acc = samples$experiment, out_type = 
##                          c("sra","submission","study","sample","experiment","run"), 
##                          sra_con = sra_con)
## write.table(conversion,file="inst/extdata/conversion.txt")


###################################################
### code chunk number 15: samples2Runs
###################################################
conversion <- read.table(system.file("extdata/conversion.txt",
                                     package="parathyroidSE",mustWork=TRUE))
samplesFull <- merge(samples, conversion)
samplesFull <- samplesFull[order(samplesFull$run),]
samplesFull <- DataFrame(lapply(samplesFull, factor))


###################################################
### code chunk number 16: addSampleColData (eval = FALSE)
###################################################
## colData(parathyroidGenesSE)$run <- sub(".*(SRR.*)_tophat_out.*","\\1",
##                                        colnames(parathyroidGenesSE))
## matchOrder <- match(colData(parathyroidGenesSE)$run, samplesFull$run)
## colData(parathyroidGenesSE) <- cbind(colData(parathyroidGenesSE),
##                                      subset(samplesFull[matchOrder,],select=-run))
## colData(parathyroidExonsSE)$run <- sub(".*(SRR.*)_tophat_out.*","\\1",
##                                        colnames(parathyroidExonsSE))
## matchOrder <- match(colData(parathyroidExonsSE)$run, samplesFull$run)
## colData(parathyroidExonsSE) <- cbind(colData(parathyroidExonsSE),
##                                      subset(samplesFull[matchOrder,],select=-run))


###################################################
### code chunk number 17: exptData (eval = FALSE)
###################################################
## exptData = new("MIAME",
##   name="Felix Haglund",
##   lab="Science for Life Laboratory Stockholm",
##   contact="Mikael Huss",
##   title="DPN and Tamoxifen treatments of parathyroid adenoma cells",
##   url="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37211",
##   abstract="Primary hyperparathyroidism (PHPT) is most frequently present in postmenopausal women. Although the involvement of estrogen has been suggested, current literature indicates that parathyroid tumors are estrogen receptor (ER) alpha negative. Objective: The aim of the study was to evaluate the expression of ERs and their putative function in parathyroid tumors. Design: A panel of 37 parathyroid tumors was analyzed for expression and promoter methylation of the ESR1 and ESR2 genes as well as expression of the ERalpha and ERbeta1/ERbeta2 proteins. Transcriptome changes in primary cultures of parathyroid adenoma cells after treatment with the selective ERbeta1 agonist diarylpropionitrile (DPN) and 4-hydroxytamoxifen were identified using next-generation RNA sequencing. Results: Immunohistochemistry revealed very low expression of ERalpha, whereas all informative tumors expressed ERbeta1 (n = 35) and ERbeta2 (n = 34). Decreased nuclear staining intensity and mosaic pattern of positive and negative nuclei of ERbeta1 were significantly associated with larger tumor size. Tumor ESR2 levels were significantly higher in female vs. male cases. In cultured cells, significantly increased numbers of genes with modified expression were detected after 48 h, compared to 24-h treatments with DPN or 4-hydroxytamoxifen, including the parathyroid-related genes CASR, VDR, JUN, CALR, and ORAI2. Bioinformatic analysis of transcriptome changes after DPN treatment revealed significant enrichment in gene sets coupled to ER activation, and a highly significant similarity to tumor cells undergoing apoptosis. Conclusions: Parathyroid tumors express ERbeta1 and ERbeta2. Transcriptional changes after ERbeta1 activation and correlation to clinical features point to a role of estrogen signaling in parathyroid function and disease.")
## pubMedIds(exptData) <- "23024189"
## exptData(parathyroidGenesSE) <- list(MIAME=exptData)
## exptData(parathyroidExonsSE) <- list(MIAME=exptData)


###################################################
### code chunk number 18: saveData (eval = FALSE)
###################################################
## save(parathyroidGenesSE,file="data/parathyroidGenesSE.RData")
## save(parathyroidExonsSE,file="data/parathyroidExonsSE.RData")


###################################################
### code chunk number 19: sessInfo
###################################################
toLatex(sessionInfo())


