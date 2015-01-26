#!/usr/bin/Rscript
library(DESeq2)
library(ggplot2)
library("RColorBrewer")
library("gplots")

args = commandArgs(TRUE)

# args = c("0.05", "3", "count.txt", "sample:extraction", "makindu:prot", "makindu:size", "chicharo:prot", "chicharo:size")
# 0.05 3 count.txt sample:extraction makindu:prot makindu:size chicharo:prot chicharo:size
if(length(args) < 4)
{
    print("countTE.R [rosette_column] [count_file] [experiment_formula] [sample_1_name ... sample_N_name]")
    exit(1)
}

FDR_level = as.numeric(args[1])

# we get the counts
count_file = args[3]
counts = read.table(args[3])
rownames(counts) = counts[,1]
counts_information = counts[, c(1:as.integer(args[2])) ]
counts_total = counts[, dim(counts)[2]]
counts = counts[, -c(1:as.integer(args[2]), dim(counts)[2]) ]

# we get the variables
variable_names = strsplit(args[4], "[:]")[[1]]
sample_names = args[5:length(args)]
variables = strsplit(sample_names[1], "[:]")[[1]]
variable_number = length(variables)

for( i in c(2:length(sample_names)))
{
    variable = strsplit(sample_names[i], "[:]")[[1]]
    variables = cbind(variables, variable)
}
variables = t(variables)
rownames(variables) = sample_names
colnames(variables) =  variable_names
variables = as.data.frame(variables)

names(counts) = sample_names
counts = as.data.frame(counts)
counts = counts[rowSums(counts)!=0,]

print(head(counts))
print(variables)

# we run DeSeq
TE = DESeqDataSetFromMatrix(countData = counts, colData = variables, design = formula(paste0("~", variable_names[1])) )

print(paste0("formula: ~", variable_names[1]))
TE = DESeq(TE, betaPrior=TRUE)

# some graphs about the quality of the analysis
res = tryCatch({
        pdf("DispEsts.pdf" , height=10,width=10)
            plotDispEsts(TE)
        dev.off()
    }, error = function(e){return(0)})

res = tryCatch({
        rld = rlogTransformation(TE, blind=T)
        data = plotPCA(rld, intgroup=variable_names, returnData=TRUE)
        percentVar = round(100 * attr(data, "percentVar"))
        ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
            geom_point(size=3) +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
            ylab(paste0("PC2: ",percentVar[2],"% variance"))
            print(plotPCA(rld, intgroup=variable_names))
        ggsave(file="PCA.pdf" , width = 10, height = 10, units = "cm")
    }, error = function(e){return(0)})

res = tryCatch({
        pdf("DispEsts.pdf" , height=10,width=10)
            plotMA(results(TE))
        dev.off()
    }, error = function(e){return(0)})

# differential analysis between every pair of variable 1
main_factor = levels(variables[,1])
i = 1
j = 1
n = length(main_factor)
number_diff = c()
for(i in c(1:n))
{
    for(j in c(1:n))
    {
        if(main_factor[i] != main_factor[j])
        {
            res = as.data.frame(results(TE, contrast=c(variable_names[1], main_factor[i], main_factor[j]), independentFiltering=F))
            res[[variable_names[1]]] = paste(main_factor[i],"vs",main_factor[j])
            res$TE = rownames(res)
            res$BH = p.adjust(res$pvalue, method="BH")
            res$signif = res$BH <= 0.1
            res = res[order(res$padj),]
            number_diff = rbind(number_diff, res)
        }
        else
        {
            tmp = as.data.frame(t(c(NA,NA,NA,NA,NA,NA)))
            names(tmp) = colnames(as.data.frame(results(TE)))
            tmp[[variable_names[1]]] = paste(main_factor[i],"vs",main_factor[j])
            tmp$TE = NA
            tmp$BH = 1
            tmp$signif = FALSE
            number_diff = rbind(number_diff, tmp)
        }
    }
}
number_diff = number_diff[!is.na(number_diff$baseMean),]

# significant FDR = 0.05
significant_TE = number_diff[number_diff$BH < FDR_level & !is.na(number_diff$log2FoldChange),]
significant_TE = significant_TE[order(significant_TE$TE),]
write.csv(number_diff, file = paste("all_TE", "~", variable_names[1], ".csv", sep=""))
write.csv(significant_TE, file = paste("significant_TE", "~", variable_names[1], ".csv", sep=""))

TE_vsd = varianceStabilizingTransformation(TE)

select = order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(assay(TE_vsd)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))


distsRL = dist(t(assay(rld)))
