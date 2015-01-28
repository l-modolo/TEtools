### R code from vignette source 'GenomicRangesHOWTOs.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: options
###################################################
options(width=72)
options("showHeadLines" = 3)
options("showTailLines" = 3)


###################################################
### code chunk number 3: load
###################################################
library(GenomicRanges)
library(Rsamtools)
library(pasillaBamSubset)
un1 <- untreated1_chr4() ## single-end reads


###################################################
### code chunk number 4: readGAlignments
###################################################
un1 <- untreated1_chr4()
gal <- readGAlignments(un1)


###################################################
### code chunk number 5: ScanBamParam
###################################################
what <- c("flag", "cigar") 
which <- GRanges("chr4", IRanges(1, 5000)) 
flag <- scanBamFlag(isMinusStrand = TRUE)
param <- ScanBamParam(which=which, what=what, flag=flag)
neg <- readGAlignments(un1, param=param)
neg


###################################################
### code chunk number 6: readGAlignmentPairs
###################################################
un3 <- untreated3_chr4()
gapairs <- readGAlignmentPairs(un3)


###################################################
### code chunk number 7: gapairs
###################################################
gapairs


###################################################
### code chunk number 8: readGAlignmentsList
###################################################
galist <- readGAlignmentsList(BamFile(un3, asMates=TRUE))


###################################################
### code chunk number 9: galist
###################################################
galist


###################################################
### code chunk number 10: non_mates
###################################################
non_mates <- galist[unlist(mcols(galist)$mates) == FALSE]
table(elementLengths(non_mates))


###################################################
### code chunk number 11: yieldSize
###################################################
bf <- BamFile(un1, yieldSize=100000)


###################################################
### code chunk number 12: readGAlignments_by_chunk
###################################################
open(bf)
cvg <- NULL
repeat {
    chunk <- readGAlignments(bf)
    if (length(chunk) == 0L)
        break
    chunk_cvg <- coverage(chunk)
    if (is.null(cvg)) {
        cvg <- chunk_cvg
    } else {
        cvg <- cvg + chunk_cvg
    }
}
close(bf)
cvg


###################################################
### code chunk number 13: load
###################################################
library(GenomicRanges)
library(Rsamtools)
library(pasillaBamSubset)
un1 <- untreated1_chr4() ## single-end records


###################################################
### code chunk number 14: count_1
###################################################
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
exbygene <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "gene")


###################################################
### code chunk number 15: count_2
###################################################
se <- summarizeOverlaps(exbygene, un1, mode="IntersectionNotEmpty")


###################################################
### code chunk number 16: count_3
###################################################
class(se)
head(table(assays(se)$counts))


###################################################
### code chunk number 17: count_4
###################################################
identical(length(exbygene), length(assays(se)$counts))


###################################################
### code chunk number 18: count_5
###################################################
rowData(se)


###################################################
### code chunk number 19: hub_1
###################################################
library(AnnotationHub)
hub <- AnnotationHub()
filters(hub) <- list(Species="Drosophila melanogaster")


###################################################
### code chunk number 20: hub_2
###################################################
length(hub)
head(names(hub))


###################################################
### code chunk number 21: hub_3
###################################################
gr <- hub$goldenpath.dm3.database.ensGene_0.0.1.RData
summary(gr)


###################################################
### code chunk number 22: hub_4
###################################################
names(metadata(gr)[[2]])
metadata(gr)[[2]]$Tags


###################################################
### code chunk number 23: hub_5
###################################################
split(gr, gr$name)


###################################################
### code chunk number 24: count_table
###################################################
library(DESeq)
deseq <- newCountDataSet(assays(se)$counts, rownames(colData(se)))
library(edgeR)
edger <- DGEList(assays(se)$counts, group=rownames(colData(se)))


###################################################
### code chunk number 25: trak_1
###################################################
trak2 <- "66008"


###################################################
### code chunk number 26: trak_2
###################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


###################################################
### code chunk number 27: trak_3
###################################################
library(GenomicFeatures)
txbygene <- transcriptsBy(txdb, by="gene")[trak2]
txbygene


###################################################
### code chunk number 28: trak_4
###################################################
tx_names <- mcols(unlist(txbygene))$tx_name
tx_names


###################################################
### code chunk number 29: trak_5
###################################################
intronsbytx <- intronsByTranscript(txdb, use.names=TRUE)[tx_names]
elementLengths(intronsbytx)


###################################################
### code chunk number 30: trak_7
###################################################
exonsbytx <- exonsBy(txdb, "tx", use.names=TRUE)[tx_names]
elementLengths(exonsbytx)


###################################################
### code chunk number 31: trak_8
###################################################
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)


###################################################
### code chunk number 32: trak_9
###################################################
intron_seqs <- extractTranscriptsFromGenome(Hsapiens, intronsbytx)
intron_seqs


###################################################
### code chunk number 33: trak_10
###################################################
exon_seqs <- extractTranscriptsFromGenome(Hsapiens, exonsbytx)
exon_seqs


###################################################
### code chunk number 34: cancer_1
###################################################
library(KEGG.db)
pathways <- toTable(KEGGPATHNAME2ID)
pathways[grepl("cancer", pathways$path_name, fixed=TRUE),] 


###################################################
### code chunk number 35: cancer_2
###################################################
library(KEGGgraph)
dest <- tempfile()
retrieveKGML("05200", "hsa", dest, "internal")


###################################################
### code chunk number 36: cancer_3
###################################################
crids <- as.character(parseKGML2DataFrame(dest)[,1])
crgenes <- unique(translateKEGGID2GeneID(crids))
head(crgenes)


###################################################
### code chunk number 37: cancer_4
###################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


###################################################
### code chunk number 38: cancer_5
###################################################
txbygene <- transcriptsBy(txdb, "gene")[crgenes] ## subset on colorectal genes
map <- relist(unlist(txbygene, use.names=FALSE)$tx_id, txbygene)
map


###################################################
### code chunk number 39: cancer_6
###################################################
cds <- cdsBy(txdb, "tx")
threeUTR <- threeUTRsByTranscript(txdb)
fiveUTR <- fiveUTRsByTranscript(txdb)


###################################################
### code chunk number 40: cancer_7
###################################################
txid <- unlist(map, use.names=FALSE)
cds <- cds[names(cds) %in% txid]
threeUTR <- threeUTR[names(threeUTR) %in% txid]
fiveUTR <- fiveUTR[names(fiveUTR) %in% txid]


###################################################
### code chunk number 41: cancer_8
###################################################
length(txid) ## all possible transcripts
length(cds)
length(threeUTR)
length(fiveUTR)


###################################################
### code chunk number 42: cancer_9
###################################################
cds


###################################################
### code chunk number 43: cancer_10
###################################################
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19


###################################################
### code chunk number 44: cancer_11
###################################################
threeUTR_seqs <- extractTranscriptsFromGenome(genome, threeUTR) 
fiveUTR_seqs <- extractTranscriptsFromGenome(genome, fiveUTR) 
cds_seqs <- extractTranscriptsFromGenome(genome, cds) 


###################################################
### code chunk number 45: cancer_12
###################################################
cds_seqs


###################################################
### code chunk number 46: cancer_13
###################################################
lst3 <- split(threeUTR_seqs, PartitioningByWidth(sum(map %in% names(threeUTR))))
lst5 <- split(fiveUTR_seqs, PartitioningByWidth(sum(map %in% names(fiveUTR))))
lstc <- split(cds_seqs, PartitioningByWidth(sum(map %in% names(cds))))
names(lst3) <- names(lst5) <- names(lstc) <- names(map)


###################################################
### code chunk number 47: cancer_14
###################################################
length(map)
table(elementLengths(map))


###################################################
### code chunk number 48: cancer_15
###################################################
table(elementLengths(lstc))
table(elementLengths(lst3))
names(lst3)[elementLengths(lst3) == 0L] ## genes with no 3' UTR data
table(elementLengths(lst5))
names(lst5)[elementLengths(lst5) == 0L] ## genes with no 5' UTR data


###################################################
### code chunk number 49: cseq_1
###################################################
library(Rsamtools)
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
param <- ScanBamParam(what=c("seq", "qual"))
gal <- readGAlignmentsFromBam(bamfile, use.names=TRUE, param=param)


###################################################
### code chunk number 50: cseq_2
###################################################
qseq <- setNames(mcols(gal)$seq, names(gal))
qual <- setNames(mcols(gal)$qual, names(gal))
qseq_on_ref <- sequenceLayer(qseq, cigar(gal),
                             from="query", to="reference")
qual_on_ref <- sequenceLayer(qual, cigar(gal),
                             from="query", to="reference")


###################################################
### code chunk number 51: cseq_3
###################################################
qseq_on_ref_by_chrom <- splitAsList(qseq_on_ref, seqnames(gal))
qual_on_ref_by_chrom <- splitAsList(qual_on_ref, seqnames(gal))
pos_by_chrom <- splitAsList(start(gal), seqnames(gal))


###################################################
### code chunk number 52: cseq_4
###################################################
gr_by_chrom <- lapply(seqlevels(gal),
  function(seqname)
  {
    qseq_on_ref2 <- qseq_on_ref_by_chrom[[seqname]]
    qual_on_ref2 <- qual_on_ref_by_chrom[[seqname]]
    pos2 <- pos_by_chrom[[seqname]]
    qseq_on_ref_per_pos <- split(qseq_on_ref2, pos2)
    qual_on_ref_per_pos <- split(qual_on_ref2, pos2)
    nread <- elementLengths(qseq_on_ref_per_pos)
    gr_mcols <- DataFrame(nread=unname(nread),
                          qseq_on_ref=unname(qseq_on_ref_per_pos),
                          qual_on_ref=unname(qual_on_ref_per_pos))
    gr <- GRanges(Rle(seqname, nrow(gr_mcols)),
                  IRanges(as.integer(names(nread)), width=1))
    mcols(gr) <- gr_mcols
    seqlevels(gr) <- seqlevels(gal)
    gr
  })


###################################################
### code chunk number 53: cseq_5
###################################################
gr <- do.call(c, gr_by_chrom)
seqinfo(gr) <- seqinfo(gal)


###################################################
### code chunk number 54: cseq_6
###################################################
gr[1:6]


###################################################
### code chunk number 55: cseq_7
###################################################
qseq_on_ref
qual_on_ref


###################################################
### code chunk number 56: cseq_8
###################################################
mcols(gr)$qseq_on_ref[[6]]


###################################################
### code chunk number 57: cseq_9
###################################################
mcols(gr)$qual_on_ref[[6]]


###################################################
### code chunk number 58: cseq_10
###################################################
qseq_on_ref <- mcols(gr)$qseq_on_ref
tmp <- unlist(qseq_on_ref, use.names=FALSE)
qseq_on_ref_id <- relist(match(tmp, tmp), qseq_on_ref)


###################################################
### code chunk number 59: cseq_11
###################################################
qseq_on_ref_id


###################################################
### code chunk number 60: cseq_12
###################################################
qseq_on_ref_id2 <- endoapply(qseq_on_ref_id,
    function(ids) ids[countMatches(ids, ids) >= 0.2 * length(ids)])


###################################################
### code chunk number 61: cseq_13
###################################################
tmp <- unlist(qseq_on_ref_id2, use.names=FALSE)
qseq_on_ref2 <- relist(unlist(qseq_on_ref, use.names=FALSE)[tmp],
                       qseq_on_ref_id2)


###################################################
### code chunk number 62: cseq_14
###################################################
split_factor <- rep.int(seqnames(gr), elementLengths(qseq_on_ref2))
qseq_on_ref2 <- unlist(qseq_on_ref2, use.names=FALSE)
qseq_on_ref2_by_chrom <- splitAsList(qseq_on_ref2, split_factor)
qseq_pos_by_chrom <- splitAsList(start(gr), split_factor)

cm_by_chrom <- lapply(names(qseq_pos_by_chrom),
    function(seqname)
        consensusMatrix(qseq_on_ref2_by_chrom[[seqname]],
                        as.prob=TRUE,
                        shift=qseq_pos_by_chrom[[seqname]]-1,
                        width=seqlengths(gr)[[seqname]]))
names(cm_by_chrom) <- names(qseq_pos_by_chrom)


###################################################
### code chunk number 63: cseq_15
###################################################
lapply(cm_by_chrom, dim)


###################################################
### code chunk number 64: cseq_16
###################################################
cs_by_chrom <- lapply(cm_by_chrom,
    function(cm) {
        ## need to "fix" 'cm' because consensusString()
        ## doesn't like consensus matrices with columns
        ## that contain only zeroes (e.g., chromosome
        ## positions with no coverage)
        idx <- colSums(cm) == 0L
        cm["+", idx] <- 1
        DNAString(consensusString(cm, ambiguityMap="N"))
    })


###################################################
### code chunk number 65: cseq_17
###################################################
cs_by_chrom


###################################################
### code chunk number 66: bim_1
###################################################
library(BSgenome.Scerevisiae.UCSC.sacCer2)
set.seed(22)
cov <- RleList(
    lapply(seqlengths(Scerevisiae),
           function(len) Rle(sample(-10:10, len, replace=TRUE))),
    compress=FALSE)
head(cov, 3)


###################################################
### code chunk number 67: bin_2
###################################################
bins1 <- tileGenome(seqinfo(Scerevisiae), tilewidth=100,
                    cut.last.tile.in.chrom=TRUE)


###################################################
### code chunk number 68: bin_3
###################################################
binnedAverage <- function(bins, numvar, mcolname)
{
    stopifnot(is(bins, "GRanges"))
    stopifnot(is(numvar, "RleList"))
    stopifnot(identical(seqlevels(bins), names(numvar)))
    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    means_list <- lapply(names(numvar),
        function(seqname) {
            views <- Views(numvar[[seqname]],
                           bins_per_chrom[[seqname]])
            viewMeans(views)
        })
    new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
    mcols(bins)[[mcolname]] <- new_mcol
    bins
}


###################################################
### code chunk number 69: bin_4
###################################################
bins1 <- binnedAverage(bins1, cov, "binned_cov")
bins1


###################################################
### code chunk number 70: SessionInfo
###################################################
sessionInfo()


