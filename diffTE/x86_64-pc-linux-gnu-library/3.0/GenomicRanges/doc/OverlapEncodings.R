### R code from vignette source 'OverlapEncodings.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: options
###################################################
options(width=100)
.precomputed_results <- system.file("doc", "precomputed_results",
                                    package="GenomicRanges", mustWork=TRUE)
.loadPrecomputed <- function(objname)
{
    filename <- paste0(objname, ".rda")
    path <- file.path(.precomputed_results, filename)
    tempenv <- new.env(parent=emptyenv())
    load(path, envir=tempenv)
    get(objname, envir=tempenv)
}
.checkIdenticalToPrecomputed <- function(obj, objname, ignore.metadata=FALSE)
{
    precomputed_obj <- .loadPrecomputed(objname)
    if (ignore.metadata)
        metadata(obj) <- metadata(precomputed_obj) <- list()
    if (!identical(obj, precomputed_obj))
        stop("'", objname, "' is not identical to precomputed version")
}


###################################################
### code chunk number 3: untreated1_chr4
###################################################
library(pasillaBamSubset)
untreated1_chr4()


###################################################
### code chunk number 4: readGAlignments
###################################################
library(GenomicRanges)
library(Rsamtools)
flag0 <- scanBamFlag(isDuplicate=FALSE, isNotPassingQualityControls=FALSE)
param0 <- ScanBamParam(flag=flag0)
U1.GAL <- readGAlignments(untreated1_chr4(), use.names=TRUE, param=param0)
head(U1.GAL)


###################################################
### code chunk number 5: U1.GAL_names_is_dup
###################################################
U1.GAL_names_is_dup <- duplicated(names(U1.GAL))
table(U1.GAL_names_is_dup)


###################################################
### code chunk number 6: U1.GAL_qnames
###################################################
U1.uqnames <- unique(names(U1.GAL))
U1.GAL_qnames <- factor(names(U1.GAL), levels=U1.uqnames)


###################################################
### code chunk number 7: U1.GAL_dup2unq
###################################################
U1.GAL_dup2unq <- match(U1.GAL_qnames, U1.GAL_qnames)


###################################################
### code chunk number 8: gaps-in-U1.GAL
###################################################
head(unique(cigar(U1.GAL)))
table(ngap(U1.GAL))


###################################################
### code chunk number 9: no-indels-in-U1.GAL
###################################################
colSums(cigarOpTable(cigar(U1.GAL)))


###################################################
### code chunk number 10: readGAlignmentPairs
###################################################
U3.galp <- readGAlignmentPairs(untreated3_chr4(), use.names=TRUE, param=param0)
head(U3.galp)


###################################################
### code chunk number 11: first-and-last-U3.galp
###################################################
head(first(U3.galp))
head(last(U3.galp))


###################################################
### code chunk number 12: isProperPair
###################################################
table(isProperPair(U3.galp))


###################################################
### code chunk number 13: keep-only-proper-pairs
###################################################
U3.GALP <- U3.galp[isProperPair(U3.galp)]


###################################################
### code chunk number 14: U3.GALP_names_is_dup
###################################################
U3.GALP_names_is_dup <- duplicated(names(U3.GALP))
table(U3.GALP_names_is_dup)


###################################################
### code chunk number 15: U3.GALP_qnames
###################################################
U3.uqnames <- unique(names(U3.GALP))
U3.GALP_qnames <- factor(names(U3.GALP), levels=U3.uqnames)


###################################################
### code chunk number 16: U3.GALP_dup2unq
###################################################
U3.GALP_dup2unq <- match(U3.GALP_qnames, U3.GALP_qnames)


###################################################
### code chunk number 17: gaps-in-U3.GALP
###################################################
head(unique(cigar(first(U3.GALP))))
head(unique(cigar(last(U3.GALP))))
table(ngap(first(U3.GALP)), ngap(last(U3.GALP)))


###################################################
### code chunk number 18: no-indels-in-U3.GALP
###################################################
colSums(cigarOpTable(cigar(first(U3.GALP))))
colSums(cigarOpTable(cigar(last(U3.GALP))))


###################################################
### code chunk number 19: U1.GAL_qseq
###################################################
param1 <- ScanBamParam(flag=flag0, what="seq", tag="NM")
U1.GAL <- readGAlignments(untreated1_chr4(), use.names=TRUE, param=param1)
U1.GAL_qseq <- mcols(U1.GAL)$seq
names(U1.GAL_qseq) <- names(U1.GAL)
head(U1.GAL_qseq)


###################################################
### code chunk number 20: U1-original-query-sequences
###################################################
U1.GAL_oqseq <- U1.GAL_qseq
U1.GAL_is_on_minus <- as.logical(strand(U1.GAL) == "-")
U1.GAL_oqseq[U1.GAL_is_on_minus] <- reverseComplement(U1.GAL_oqseq[U1.GAL_is_on_minus])
head(U1.GAL_oqseq)


###################################################
### code chunk number 21: same-name-implies-same-seq-in-U1-oqseq
###################################################
stopifnot(all(U1.GAL_oqseq == U1.GAL_oqseq[U1.GAL_dup2unq]))


###################################################
### code chunk number 22: U1.oqseq
###################################################
U1.oqseq <- U1.GAL_oqseq[!U1.GAL_names_is_dup]


###################################################
### code chunk number 23: Dmelanogaster
###################################################
library(BSgenome.Dmelanogaster.UCSC.dm3)
Dmelanogaster


###################################################
### code chunk number 24: U1.grl
###################################################
U1.grl <- grglist(U1.GAL, order.as.in.query=TRUE)


###################################################
### code chunk number 25: U1-reference-query-sequences
###################################################
library(GenomicFeatures)
U1.GAL_rqseq <- extractTranscriptsFromGenome(Dmelanogaster, U1.grl)
head(U1.GAL_rqseq)


###################################################
### code chunk number 26: U1-oqseq-vs-U1-rqseq
###################################################
U1.GAL_nedit500 <- sapply(1:500, function(i) neditAt(U1.GAL_oqseq[[i]], U1.GAL_rqseq[[i]]))
table(U1.GAL_nedit500)


###################################################
### code chunk number 27: U1.GAL-nedit-vs-NM
###################################################
U1.GAL_NM <- mcols(U1.GAL)$NM
stopifnot(all(U1.GAL_NM[1:500] == U1.GAL_nedit500))


###################################################
### code chunk number 28: up-to-6-mismatches
###################################################
table(U1.GAL_NM)


###################################################
### code chunk number 29: reload-U3.GALP
###################################################
flag2 <- scanBamFlag(isDuplicate=FALSE,
                     isNotPassingQualityControls=FALSE,
                     isProperPair=TRUE)
param2 <- ScanBamParam(flag=flag2, what="seq", tag="NM")
U3.GALP <- readGAlignmentPairs(untreated3_chr4(), use.names=TRUE, param=param2)


###################################################
### code chunk number 30: U3.GALP-qseq1-and-qseq2
###################################################
U3.GALP_qseq1 <- mcols(first(U3.GALP))$seq
U3.GALP_qseq2 <- mcols(last(U3.GALP))$seq
names(U3.GALP_qseq1) <- names(U3.GALP_qseq2) <- names(U3.GALP)
head(U3.GALP_qseq1)
head(U3.GALP_qseq2)


###################################################
### code chunk number 31: U3-original-query-sequences
###################################################
U3.GALP_oqseq1 <- U3.GALP_qseq1
U3.GALP_first_is_on_minus <- as.logical(strand(first(U3.GALP)) == "-")
U3.GALP_oqseq1[U3.GALP_first_is_on_minus] <- reverseComplement(U3.GALP_oqseq1[U3.GALP_first_is_on_minus])

U3.GALP_oqseq2 <- U3.GALP_qseq2
U3.GALP_last_is_on_minus <- as.logical(strand(last(U3.GALP)) == "-")
U3.GALP_oqseq2[U3.GALP_last_is_on_minus] <- reverseComplement(U3.GALP_oqseq2[U3.GALP_last_is_on_minus])


###################################################
### code chunk number 32: same-name-implies-same-seq-pair-in-U3-oqseq
###################################################
stopifnot(all(U3.GALP_oqseq1 == U3.GALP_oqseq1[U3.GALP_dup2unq]))
stopifnot(all(U3.GALP_oqseq2 == U3.GALP_oqseq2[U3.GALP_dup2unq]))


###################################################
### code chunk number 33: U3.oqseq
###################################################
U3.oqseq1 <- U3.GALP_oqseq1[!U3.GALP_names_is_dup]
U3.oqseq2 <- U3.GALP_oqseq2[!U3.GALP_names_is_dup]


###################################################
### code chunk number 34: U3.grl_first-and-U3.grl_last
###################################################
U3.grl_first <- grglist(first(U3.GALP), order.as.in.query=TRUE)
U3.grl_last <- grglist(last(U3.GALP, invert.strand=TRUE), order.as.in.query=TRUE)


###################################################
### code chunk number 35: U3-reference-query-sequences
###################################################
U3.GALP_rqseq1 <- extractTranscriptsFromGenome(Dmelanogaster, U3.grl_first)
U3.GALP_rqseq2 <- extractTranscriptsFromGenome(Dmelanogaster, U3.grl_last)


###################################################
### code chunk number 36: U3-oqseq-vs-U3-rqseq
###################################################
U3.GALP_first_nedit500 <- sapply(1:500, function(i)
    neditAt(U3.GALP_oqseq1[[i]], U3.GALP_rqseq1[[i]])
)
table(U3.GALP_first_nedit500)

U3.GALP_last_nedit500 <- sapply(1:500, function(i)
    neditAt(reverseComplement(U3.GALP_oqseq2[[i]]), U3.GALP_rqseq2[[i]])
)
table(U3.GALP_last_nedit500)


###################################################
### code chunk number 37: U3.GALP-nedit-vs-NM
###################################################
U3.GALP_first_NM <- mcols(first(U3.GALP))$NM
stopifnot(all(U3.GALP_first_NM[1:500] == U3.GALP_first_nedit500))

U3.GALP_last_NM <- mcols(last(U3.GALP))$NM
stopifnot(all(U3.GALP_last_NM[1:500] == U3.GALP_last_nedit500))


###################################################
### code chunk number 38: up-to-2-mismatches-per-end
###################################################
table(U3.GALP_first_NM, U3.GALP_last_NM)


###################################################
### code chunk number 39: txdb
###################################################
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
TxDb.Dmelanogaster.UCSC.dm3.ensGene
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene


###################################################
### code chunk number 40: exbytx
###################################################
exbytx <- exonsBy(txdb, by="tx", use.names=TRUE)
length(exbytx)  # nb of transcripts


###################################################
### code chunk number 41: CHECK_exbytx
###################################################
.checkIdenticalToPrecomputed(exbytx, "exbytx", ignore.metadata=TRUE)


###################################################
### code chunk number 42: check-for-trans-splicing-in-exbytx
###################################################
table(elementLengths(runLength(seqnames(exbytx))))
table(elementLengths(runLength(strand(exbytx))))


###################################################
### code chunk number 43: exbytx_strand
###################################################
exbytx_strand <- unlist(runValue(strand(exbytx)), use.names=FALSE)


###################################################
### code chunk number 44: exbytx2gene
###################################################
tx <- transcripts(txdb, columns=c("tx_name", "gene_id"))
head(tx)
df <- mcols(tx)
exbytx2gene <- as.character(df$gene_id)
exbytx2gene <- factor(exbytx2gene, levels=unique(exbytx2gene))
names(exbytx2gene) <- df$tx_name
exbytx2gene <- exbytx2gene[names(exbytx)]
head(exbytx2gene)
nlevels(exbytx2gene)  # nb of genes


###################################################
### code chunk number 45: U1.OV00
###################################################
U1.OV00 <- findOverlaps(U1.GAL, exbytx, ignore.strand=TRUE)


###################################################
### code chunk number 46: length-of-U1.OV00
###################################################
length(U1.OV00)


###################################################
### code chunk number 47: nhitPerQuery-and-nhitPerSubject
###################################################
nhitPerQuery <- function(x) tabulate(queryHits(x), nbins=queryLength(x))
nhitPerSubject <- function(x) tabulate(subjectHits(x), nbins=subjectLength(x))


###################################################
### code chunk number 48: U1.GAL_ntx
###################################################
U1.GAL_ntx <- nhitPerQuery(U1.OV00)
mcols(U1.GAL)$ntx <- U1.GAL_ntx
head(U1.GAL)
table(U1.GAL_ntx)
mean(U1.GAL_ntx >= 1)


###################################################
### code chunk number 49: U1.GAL_ntx_again (eval = FALSE)
###################################################
## U1.GAL_ntx_again <- countOverlaps(U1.GAL, exbytx, ignore.strand=TRUE)
## stopifnot(identical(unname(U1.GAL_ntx_again), U1.GAL_ntx))


###################################################
### code chunk number 50: U1.uqnames_ntx
###################################################
U1.OV10 <- remapHits(U1.OV00, query.map=U1.GAL_qnames)
U1.uqnames_ntx <- nhitPerQuery(U1.OV10)
names(U1.uqnames_ntx) <- U1.uqnames
table(U1.uqnames_ntx)
mean(U1.uqnames_ntx >= 1)


###################################################
### code chunk number 51: U1.exbytx_nOV10
###################################################
U1.exbytx_nOV10 <- nhitPerSubject(U1.OV10)
names(U1.exbytx_nOV10) <- names(exbytx)
mean(U1.exbytx_nOV10 >= 50)


###################################################
### code chunk number 52: top-10-transcripts-based-on-U1.exbytx_nOV10
###################################################
head(sort(U1.exbytx_nOV10, decreasing=TRUE), n=10) 


###################################################
### code chunk number 53: U3.OV00
###################################################
U3.OV00 <- findOverlaps(U3.GALP, exbytx, ignore.strand=TRUE)


###################################################
### code chunk number 54: length-of-U3.OV00
###################################################
length(U3.OV00)


###################################################
### code chunk number 55: U3.GALP_ntx
###################################################
U3.GALP_ntx <- nhitPerQuery(U3.OV00)
mcols(U3.GALP)$ntx <- U3.GALP_ntx
head(U3.GALP)
table(U3.GALP_ntx)
mean(U3.GALP_ntx >= 1)


###################################################
### code chunk number 56: U3.GALP_ntx_again (eval = FALSE)
###################################################
## U3.GALP_ntx_again <- countOverlaps(U3.GALP, exbytx, ignore.strand=TRUE)
## stopifnot(identical(unname(U3.GALP_ntx_again), U3.GALP_ntx))


###################################################
### code chunk number 57: U3.uqnames_ntx
###################################################
U3.OV10 <- remapHits(U3.OV00, query.map=U3.GALP_qnames)
U3.uqnames_ntx <- nhitPerQuery(U3.OV10)
names(U3.uqnames_ntx) <- U3.uqnames
table(U3.uqnames_ntx)
mean(U3.uqnames_ntx >= 1)


###################################################
### code chunk number 58: U3.exbytx_nOV10
###################################################
U3.exbytx_nOV10 <- nhitPerSubject(U3.OV10)
names(U3.exbytx_nOV10) <- names(exbytx)
mean(U3.exbytx_nOV10 >= 50)


###################################################
### code chunk number 59: top-10-transcripts-based-on-U3.exbytx_nOV10
###################################################
head(sort(U3.exbytx_nOV10, decreasing=TRUE), n=10)


###################################################
### code chunk number 60: U1.grlf
###################################################
U1.grlf <- flipQuery(U1.grl)  # flipped


###################################################
### code chunk number 61: U1.ovencAB
###################################################
U1.ovencA <- encodeOverlaps(U1.grl, exbytx, hits=U1.OV00)
U1.ovencB <- encodeOverlaps(U1.grlf, exbytx, hits=U1.OV00)


###################################################
### code chunk number 62: U1.ovenc
###################################################
U1.grl_strand <- unlist(runValue(strand(U1.grl)), use.names=FALSE)
U1.ovenc <- selectEncodingWithCompatibleStrand(U1.ovencA, U1.ovencB,
                                               U1.grl_strand, exbytx_strand,
                                               hits=U1.OV00)
U1.ovenc


###################################################
### code chunk number 63: U1.ovenc_again
###################################################
U1.ovenc_again <- encodeOverlaps(U1.grl, exbytx, hits=U1.OV00, flip.query.if.wrong.strand=TRUE)
stopifnot(identical(U1.ovenc_again, U1.ovenc))


###################################################
### code chunk number 64: U1.ovenc_table
###################################################
U1.unique_encodings <- levels(U1.ovenc)
length(U1.unique_encodings)
head(U1.unique_encodings)
U1.ovenc_table <- table(encoding(U1.ovenc))
tail(sort(U1.ovenc_table))


###################################################
### code chunk number 65: U3.ovenc
###################################################
U3.grl <- grglist(U3.GALP, order.as.in.query=TRUE)
U3.ovenc <- encodeOverlaps(U3.grl, exbytx, hits=U3.OV00, flip.query.if.wrong.strand=TRUE)
U3.ovenc


###################################################
### code chunk number 66: U3.ovenc_table
###################################################
U3.unique_encodings <- levels(U3.ovenc)
length(U3.unique_encodings)
head(U3.unique_encodings)
U3.ovenc_table <- table(encoding(U3.ovenc))
tail(sort(U3.ovenc_table))


###################################################
### code chunk number 67: U1-unique-compatible-encodings
###################################################
sort(U1.ovenc_table[isCompatibleWithSplicing(U1.unique_encodings)])


###################################################
### code chunk number 68: U1.OV00_is_comp
###################################################
U1.OV00_is_comp <- isCompatibleWithSplicing(U1.ovenc)
table(U1.OV00_is_comp)  # 531797 "compatible" overlaps


###################################################
### code chunk number 69: U1.compOV00
###################################################
U1.compOV00 <- U1.OV00[U1.OV00_is_comp]


###################################################
### code chunk number 70: U1.compOV00_again (eval = FALSE)
###################################################
## U1.compOV00_again <- findCompatibleOverlaps(U1.GAL, exbytx)
## stopifnot(identical(U1.compOV00_again, U1.compOV00))


###################################################
### code chunk number 71: U1.GAL_ncomptx
###################################################
U1.GAL_ncomptx <- nhitPerQuery(U1.compOV00)
mcols(U1.GAL)$ncomptx <- U1.GAL_ncomptx
head(U1.GAL)
table(U1.GAL_ncomptx)
mean(U1.GAL_ncomptx >= 1)


###################################################
### code chunk number 72: U1.GAL_ncomptx_again (eval = FALSE)
###################################################
## U1.GAL_ncomptx_again <- countCompatibleOverlaps(U1.GAL, exbytx)
## stopifnot(identical(U1.GAL_ncomptx_again, U1.GAL_ncomptx))


###################################################
### code chunk number 73: U1.uqnames_ncomptx
###################################################
U1.compOV10 <- remapHits(U1.compOV00, query.map=U1.GAL_qnames)
U1.uqnames_ncomptx <- nhitPerQuery(U1.compOV10)
names(U1.uqnames_ncomptx) <- U1.uqnames
table(U1.uqnames_ncomptx)
mean(U1.uqnames_ncomptx >= 1)


###################################################
### code chunk number 74: U1.exbytx_ncompOV10
###################################################
U1.exbytx_ncompOV10 <- nhitPerSubject(U1.compOV10)
names(U1.exbytx_ncompOV10) <- names(exbytx)
mean(U1.exbytx_ncompOV10 >= 50)


###################################################
### code chunk number 75: top-10-transcripts-based-on-U1.exbytx_ncompOV10
###################################################
head(sort(U1.exbytx_ncompOV10, decreasing=TRUE), n=10)


###################################################
### code chunk number 76: U3-unique-compatible-encodings
###################################################
sort(U3.ovenc_table[isCompatibleWithSplicing(U3.unique_encodings)])


###################################################
### code chunk number 77: U3.OV00_is_comp
###################################################
U3.OV00_is_comp <- isCompatibleWithSplicing(U3.ovenc)
table(U3.OV00_is_comp)  # 106835 "compatible" paired-end overlaps


###################################################
### code chunk number 78: U3.compOV00
###################################################
U3.compOV00 <- U3.OV00[U3.OV00_is_comp]


###################################################
### code chunk number 79: U3.compOV00_again (eval = FALSE)
###################################################
## U3.compOV00_again <- findCompatibleOverlaps(U3.GALP, exbytx)
## stopifnot(identical(U3.compOV00_again, U3.compOV00))


###################################################
### code chunk number 80: U3.GALP_ncomptx
###################################################
U3.GALP_ncomptx <- nhitPerQuery(U3.compOV00)
mcols(U3.GALP)$ncomptx <- U3.GALP_ncomptx
head(U3.GALP)
table(U3.GALP_ncomptx)
mean(U3.GALP_ncomptx >= 1)


###################################################
### code chunk number 81: U3.GALP_ncomptx_again (eval = FALSE)
###################################################
## U3.GALP_ncomptx_again <- countCompatibleOverlaps(U3.GALP, exbytx)
## stopifnot(identical(U3.GALP_ncomptx_again, U3.GALP_ncomptx))


###################################################
### code chunk number 82: U3.uqnames_ncomptx
###################################################
U3.compOV10 <- remapHits(U3.compOV00, query.map=U3.GALP_qnames)
U3.uqnames_ncomptx <- nhitPerQuery(U3.compOV10)
names(U3.uqnames_ncomptx) <- U3.uqnames
table(U3.uqnames_ncomptx)
mean(U3.uqnames_ncomptx >= 1)


###################################################
### code chunk number 83: U3.exbytx_ncompOV10
###################################################
U3.exbytx_ncompOV10 <- nhitPerSubject(U3.compOV10)
names(U3.exbytx_ncompOV10) <- names(exbytx)
mean(U3.exbytx_ncompOV10 >= 50)


###################################################
### code chunk number 84: top-10-transcripts-based-on-U3.exbytx_ncompOV10
###################################################
head(sort(U3.exbytx_ncompOV10, decreasing=TRUE), n=10)


###################################################
### code chunk number 85: U1.OV00_qstart
###################################################
U1.OV00_qstart <- extractQueryStartInTranscript(U1.grl, exbytx,
                                                hits=U1.OV00, ovenc=U1.ovenc)
head(subset(U1.OV00_qstart, U1.OV00_is_comp))


###################################################
### code chunk number 86: txseq
###################################################
txseq <- extractTranscriptsFromGenome(Dmelanogaster, exbytx)


###################################################
### code chunk number 87: U1.OV00_rqseq-vs-U1.OV00_txseq
###################################################
U1.OV00_rqseq <- U1.GAL_rqseq[queryHits(U1.OV00)]
U1.OV00_rqseq[flippedQuery(U1.ovenc)] <- reverseComplement(U1.OV00_rqseq[flippedQuery(U1.ovenc)])
U1.OV00_txseq <- txseq[subjectHits(U1.OV00)]
stopifnot(all(
    U1.OV00_rqseq[U1.OV00_is_comp] ==
        narrow(U1.OV00_txseq[U1.OV00_is_comp],
               start=U1.OV00_qstart$startInTranscript[U1.OV00_is_comp],
               width=width(U1.OV00_rqseq)[U1.OV00_is_comp])
))


###################################################
### code chunk number 88: U3.OV00_Lqstart
###################################################
U3.OV00_Lqstart <- extractQueryStartInTranscript(U3.grl, exbytx,
                                                 hits=U3.OV00, ovenc=U3.ovenc)
head(subset(U3.OV00_Lqstart, U3.OV00_is_comp))


###################################################
### code chunk number 89: U3.OV00_Rqstart
###################################################
U3.OV00_Rqstart <- extractQueryStartInTranscript(U3.grl, exbytx,
                                                 hits=U3.OV00, ovenc=U3.ovenc,
                                                 for.query.right.end=TRUE)
head(subset(U3.OV00_Rqstart, U3.OV00_is_comp))


###################################################
### code chunk number 90: U3.OV00_Lrqseq_and_Rrqseq
###################################################
U3.OV00_Lrqseq <- U3.GALP_rqseq1[queryHits(U3.OV00)]
U3.OV00_Rrqseq <- U3.GALP_rqseq2[queryHits(U3.OV00)]


###################################################
### code chunk number 91: U3.OV00_Lrqseq_and_Rrqseq
###################################################
flip_idx <- which(flippedQuery(U3.ovenc))
tmp <- U3.OV00_Lrqseq[flip_idx]
U3.OV00_Lrqseq[flip_idx] <- reverseComplement(U3.OV00_Rrqseq[flip_idx])
U3.OV00_Rrqseq[flip_idx] <- reverseComplement(tmp)


###################################################
### code chunk number 92: U3.OV00_txseq
###################################################
U3.OV00_txseq <- txseq[subjectHits(U3.OV00)]


###################################################
### code chunk number 93: U3.OV00_Lrqseq-vs-U3.OV00_txseq
###################################################
stopifnot(all(
    U3.OV00_Lrqseq[U3.OV00_is_comp] ==
        narrow(U3.OV00_txseq[U3.OV00_is_comp],
               start=U3.OV00_Lqstart$startInTranscript[U3.OV00_is_comp],
               width=width(U3.OV00_Lrqseq)[U3.OV00_is_comp])
))


###################################################
### code chunk number 94: U3.OV00_Rrqseq-vs-U3.OV00_txseq
###################################################
stopifnot(all(
    U3.OV00_Rrqseq[U3.OV00_is_comp] ==
        narrow(U3.OV00_txseq[U3.OV00_is_comp],
               start=U3.OV00_Rqstart$startInTranscript[U3.OV00_is_comp],
               width=width(U3.OV00_Rrqseq)[U3.OV00_is_comp])
))


###################################################
### code chunk number 95: findSequenceHits
###################################################
### A wrapper to vwhichPDict() that supports IUPAC ambiguity codes in 'qseq'
### and 'txseq', and treats them as such.
findSequenceHits <- function(qseq, txseq, which.txseq=NULL, max.mismatch=0)
{
    .asHits <- function(x, pattern_length)
    {
        query_hits <- unlist(x)
        if (is.null(query_hits))
            query_hits <- integer(0)
        subject_hits <- rep.int(seq_len(length(x)), elementLengths(x))
        new("Hits", queryHits=query_hits, subjectHits=subject_hits,
                    queryLength=pattern_length, subjectLength=length(x))
    }

    .isHitInTranscriptBounds <- function(hits, qseq, txseq)
    {
        sapply(seq_len(length(hits)),
            function(i) {
                pattern <- qseq[[queryHits(hits)[i]]]
                subject <- txseq[[subjectHits(hits)[i]]]
                v <- matchPattern(pattern, subject,
                                  max.mismatch=max.mismatch, fixed=FALSE)
                any(1L <= start(v) & end(v) <= length(subject))
            })
    }

    if (!is.null(which.txseq)) {
        txseq0 <- txseq
        txseq <- txseq[which.txseq]
    }

    names(qseq) <- NULL
    other <- alphabetFrequency(qseq, baseOnly=TRUE)[ , "other"]
    is_clean <- other == 0L  # "clean" means "no IUPAC ambiguity code"

    ## Find hits for "clean" original queries.
    qseq0 <- qseq[is_clean]
    pdict0 <- PDict(qseq0, max.mismatch=max.mismatch)
    m0 <- vwhichPDict(pdict0, txseq,
                      max.mismatch=max.mismatch, fixed="pattern")
    hits0 <- .asHits(m0, length(qseq0))
    hits0@queryLength <- length(qseq)
    hits0@queryHits <- which(is_clean)[hits0@queryHits]

    ## Find hits for non "clean" original queries.
    qseq1 <- qseq[!is_clean]
    m1 <- vwhichPDict(qseq1, txseq,
                      max.mismatch=max.mismatch, fixed=FALSE)
    hits1 <- .asHits(m1, length(qseq1))
    hits1@queryLength <- length(qseq)
    hits1@queryHits <- which(!is_clean)[hits1@queryHits]

    ## Combine the hits.
    query_hits <- c(queryHits(hits0), queryHits(hits1))
    subject_hits <- c(subjectHits(hits0), subjectHits(hits1))

    if (!is.null(which.txseq)) {
        ## Remap the hits.
        txseq <- txseq0
        subject_hits <- which.txseq[subject_hits]
        hits0@subjectLength <- length(txseq)
    }

    ## Order the hits.
    oo <- IRanges:::orderIntegerPairs(query_hits, subject_hits)
    hits0@queryHits <- query_hits[oo]
    hits0@subjectHits <- subject_hits[oo]

    if (max.mismatch != 0L) {
        ## Keep only "in bounds" hits.
        is_in_bounds <- .isHitInTranscriptBounds(hits0, qseq, txseq)
        hits0 <- hits0[is_in_bounds]
    }
    hits0
}


###################################################
### code chunk number 96: which.txseq (eval = FALSE)
###################################################
## chr4tx <- transcripts(txdb, vals=list(tx_chrom="chr4"))
## chr4txnames <- mcols(chr4tx)$tx_name
## which.txseq <- match(chr4txnames, names(txseq))


###################################################
### code chunk number 97: U1.sbcompHITS (eval = FALSE)
###################################################
## U1.sbcompHITSa <- findSequenceHits(U1.oqseq, txseq,
##                                    which.txseq=which.txseq, max.mismatch=6)
## U1.sbcompHITSb <- findSequenceHits(reverseComplement(U1.oqseq), txseq,
##                                    which.txseq=which.txseq, max.mismatch=6)
## U1.sbcompHITS <- union(U1.sbcompHITSa, U1.sbcompHITSb)


###################################################
### code chunk number 98: LOAD_U1.sbcompHITS
###################################################
U1.sbcompHITSa <- .loadPrecomputed("U1.sbcompHITSa")
U1.sbcompHITSb <- .loadPrecomputed("U1.sbcompHITSb")
U1.sbcompHITS <- union(U1.sbcompHITSa, U1.sbcompHITSb)


###################################################
### code chunk number 99: U1.uqnames_nsbcomptx
###################################################
U1.uqnames_nsbcomptx <- nhitPerQuery(U1.sbcompHITS)
names(U1.uqnames_nsbcomptx) <- U1.uqnames
table(U1.uqnames_nsbcomptx)
mean(U1.uqnames_nsbcomptx >= 1)


###################################################
### code chunk number 100: U1.exbytx_nsbcompHITS
###################################################
U1.exbytx_nsbcompHITS <- nhitPerSubject(U1.sbcompHITS)
names(U1.exbytx_nsbcompHITS) <- names(exbytx)
mean(U1.exbytx_nsbcompHITS >= 50)


###################################################
### code chunk number 101: top-10-transcripts-based-on-U1.exbytx_nsbcompHITS
###################################################
head(sort(U1.exbytx_nsbcompHITS, decreasing=TRUE), n=10)


###################################################
### code chunk number 102: encoding-based-compatible-implies-string-based-compatible
###################################################
stopifnot(length(setdiff(U1.compOV10, U1.sbcompHITS)) == 0)


###################################################
### code chunk number 103: string-based-compatible-does-NOT-imply-encoding-based-compatible
###################################################
length(setdiff(U1.sbcompHITS, U1.compOV10))


###################################################
### code chunk number 104: U1-unique-almost-compatible-encodings
###################################################
sort(U1.ovenc_table[isCompatibleWithSkippedExons(U1.unique_encodings)])


###################################################
### code chunk number 105: U1.OV00_is_acomp
###################################################
U1.OV00_is_acomp <- isCompatibleWithSkippedExons(U1.ovenc)
table(U1.OV00_is_acomp)  # 1202 "almost compatible" overlaps


###################################################
### code chunk number 106: U1.acompOV00
###################################################
U1.acompOV00 <- U1.OV00[U1.OV00_is_acomp]


###################################################
### code chunk number 107: U1.GAL_nacomptx
###################################################
U1.GAL_nacomptx <- nhitPerQuery(U1.acompOV00)
mcols(U1.GAL)$nacomptx <- U1.GAL_nacomptx
head(U1.GAL)
table(U1.GAL_nacomptx)
mean(U1.GAL_nacomptx >= 1)


###################################################
### code chunk number 108: U1.exbytx_nacompOV00
###################################################
U1.exbytx_nacompOV00 <- nhitPerSubject(U1.acompOV00)
names(U1.exbytx_nacompOV00) <- names(exbytx)
table(U1.exbytx_nacompOV00)
mean(U1.exbytx_nacompOV00 >= 50)


###################################################
### code chunk number 109: U1.OV00_qstart
###################################################
head(subset(U1.OV00_qstart, U1.OV00_is_acomp))


###################################################
### code chunk number 110: U3-unique-almost-compatible-encodings
###################################################
sort(U3.ovenc_table[isCompatibleWithSkippedExons(U3.unique_encodings)])


###################################################
### code chunk number 111: U3.OV00_is_acomp
###################################################
U3.OV00_is_acomp <- isCompatibleWithSkippedExons(U3.ovenc)
table(U3.OV00_is_acomp)  # 141 "almost compatible" paired-end overlaps


###################################################
### code chunk number 112: U3.acompOV00
###################################################
U3.acompOV00 <- U3.OV00[U3.OV00_is_acomp]


###################################################
### code chunk number 113: U3.GALP_nacomptx
###################################################
U3.GALP_nacomptx <- nhitPerQuery(U3.acompOV00)
mcols(U3.GALP)$nacomptx <- U3.GALP_nacomptx
head(U3.GALP)
table(U3.GALP_nacomptx)
mean(U3.GALP_nacomptx >= 1)


###################################################
### code chunk number 114: U3.exbytx_nacompOV00
###################################################
U3.exbytx_nacompOV00 <- nhitPerSubject(U3.acompOV00)
names(U3.exbytx_nacompOV00) <- names(exbytx)
table(U3.exbytx_nacompOV00)
mean(U3.exbytx_nacompOV00 >= 50)


###################################################
### code chunk number 115: U3.OV00_Lqstart-and-U3.OV00_Rqstart
###################################################
head(subset(U3.OV00_Lqstart, U3.OV00_is_acomp))
head(subset(U3.OV00_Rqstart, U3.OV00_is_acomp))


###################################################
### code chunk number 116: U1.GAL_is_nsj
###################################################
U1.GAL_is_nsj <- U1.GAL_nacomptx != 0L & U1.GAL_ncomptx == 0L
head(which(U1.GAL_is_nsj))


###################################################
### code chunk number 117: U1.OV00_is_nsj
###################################################
U1.OV00_is_nsj <- queryHits(U1.OV00) %in% which(U1.GAL_is_nsj)


###################################################
### code chunk number 118: narrow-U1.OV00_is_nsj
###################################################
U1.OV00_is_nsj <- U1.OV00_is_nsj & U1.OV00_is_acomp
U1.nsjOV00 <- U1.OV00[U1.OV00_is_nsj]


###################################################
### code chunk number 119: U1.nsjOV00_skippedex
###################################################
U1.nsjOV00_skippedex <- extractSkippedExonRanks(U1.ovenc)[U1.OV00_is_nsj]
names(U1.nsjOV00_skippedex) <- queryHits(U1.nsjOV00)
table(elementLengths(U1.nsjOV00_skippedex))


###################################################
### code chunk number 120: U1.exbytx_skippedex
###################################################
f <- factor(names(exbytx)[subjectHits(U1.nsjOV00)], levels=names(exbytx))
U1.exbytx_skippedex <- split(U1.nsjOV00_skippedex, f)


###################################################
### code chunk number 121: names-of-U1.exbytx_skippedex
###################################################
head(names(U1.exbytx_skippedex))  # transcript names


###################################################
### code chunk number 122: FBtr0089124-skipped-exons
###################################################
U1.exbytx_skippedex$FBtr0089124


###################################################
### code chunk number 123: FBtr0089147-skipped-exons
###################################################
U1.exbytx_skippedex$FBtr0089147


###################################################
### code chunk number 124: sessionInfo
###################################################
sessionInfo()


