#!/usr/bin/Rscript

# sudo apt-get install libxml2-dev
# yum install libxml2-devel glibc
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2", dep=T)
# biocLite("gplots", dep=T)
# biocLite("ggplot2", dep=T)
# biocLite("RColorBrewer", dep=T)

args = commandArgs(trailingOnly=FALSE)
scriptPath <- dirname(sub("--file=","",args[grep("--file",args)]))
scriptPath = paste0(scriptPath,"/diffTE")
print(scriptPath)
.libPaths(c(.libPaths(), scriptPath))

args = commandArgs(trailingOnly=TRUE)
for(i in 2:length(args))
{
    eval(parse(text=sub("--", "", args[i])))
}

help = FALSE
if(!exists("FDR_level")){help = TRUE; print("FDR_level")}
if(!exists("count_column")){help = TRUE; print("count_column")}
if(!exists("count_file")){help = TRUE; print("count_file")}
if(!exists("experiment_formula")){help = TRUE; print("experiment_formula")}
if(!exists("sample_names")){help = TRUE; print("sample_names")}

if(help==TRUE)
{
    print("diffTE.R --args --FDR_level=0.05 --count_column=2 --count_file=\\\"count.txt\\\" --experiment_formula=\\\"population:type\\\" --sample_names=\\\"population1:type1,population1:type2,population2:type1,population2:type2\\\"")
    quit("no")
}
if(exists("version"))
{
    print("1.0.0")
}
count_column = count_column-1
print(FDR_level)
print(count_column)
print(count_file)
print(experiment_formula)
print(sample_names)

# we get the counts
counts = read.table(count_file)
rownames(counts) = counts[,1]
counts_information = counts[, c(1:count_column) ]
counts_total = counts[, dim(counts)[2]]
counts = counts[, -c(1:count_column, dim(counts)[2]) ]

# we get the variables
variable_names = strsplit(experiment_formula, "[:]")[[1]]
sample_names = strsplit(sample_names,",")[[1]]
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
suppressMessages(require(genefilter, quietly = TRUE))
suppressMessages(require(DESeq2, quietly = TRUE))
suppressMessages(require(gplots, quietly = TRUE))
suppressMessages(require(ggplot2, quietly = TRUE))
suppressMessages(require(RColorBrewer, quietly = TRUE))


TE = DESeqDataSetFromMatrix(countData = counts, 
                            colData = variables, 
                            design = formula(paste0("~", variable_names[1]))
                            )

print(paste0("formula: ~", variable_names[1]))
TE = DESeq(TE, betaPrior=TRUE)

# some graphs about the quality of the analysis
pdf("DispEsts.pdf" , height=10,width=10)
    plotDispEsts(TE)
x = dev.off()

ntop = 500
rld = rlogTransformation(TE, blind=T)
rv = rowVars(assay(rld))
select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca = prcomp(t(assay(rld)[select, ]))

if(dim(variables)[2] == 2)
{
    ggplot(data=as.data.frame(pca$x), aes(PC1, PC2, color=variables[,1], shape=variables[,2] )) + 
            geom_point(size = 6) + 
            xlab(paste0("PC1: ",100*summary(pca)[6]$importance[2,][1],"% of variance")) + 
            ylab(paste0("PC1: ",100*summary(pca)[6]$importance[2,][2],"% of variance")) + 
            theme_bw() + 
            guides(color=guide_legend(title=variable_names[1]), shape=guide_legend(title=variable_names[2]))
    ggsave(file="PCA.pdf", width=20, height=20, units="cm", dpi=1200)
} else {
    fac = factor(apply(as.data.frame(colData(rld)[, variable_names, drop = FALSE]),
        1, paste, collapse = " : "))
    ggplot(data=as.data.frame(pca$x), aes(PC1, PC2, color=fac)) + 
            geom_point(size = 6) + 
            xlab(paste0("PC1: ",100*summary(pca)[6]$importance[2,][1],"% of variance")) + 
            ylab(paste0("PC1: ",100*summary(pca)[6]$importance[2,][2],"% of variance")) + 
            theme_bw() + 
            guides(color=guide_legend(title="factors"))
    ggsave(file="PCA.pdf", width=20, height=20, units="cm", dpi=1200)
}

pdf("MA.pdf" , height=10,width=10)
    plotMA(results(TE))
x = dev.off()

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
            res = as.data.frame(results(TE, 
                                        contrast=c(variable_names[1], 
                                                    main_factor[i], 
                                                    main_factor[j]), 
                                        independentFiltering=F
                                        )
                                )
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
write.csv(number_diff, file = paste0("all_TE", "~", variable_names[1], ".csv"))
write.csv(significant_TE, file = paste0("significant_TE", "~", variable_names[1], ".csv"))

TE_vsd = varianceStabilizingTransformation(TE)
TE_row = order(rowMeans(counts(TE,normalized=TRUE)),decreasing=TRUE)
old_i = 1
for(i in seq(from=30, to=length(TE_row), by = 30))
{
    select = order(rownames(TE),decreasing=FALSE)[old_i:i]
    hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    pdf(paste0("heatmap_", old_i, "-", i, ".pdf") , height=10,width=10)
    heatmap.2(assay(TE_vsd)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
    old_i = i
}

# volcanoplot
ggplot(number_diff, aes(x=log2FoldChange, y=-log2(BH), colour = BH)) +
        geom_point(size = 1) +
        facet_wrap(as.formula(paste0("~", variable_names[1]))) +
         xlab("log2 foldchange") +
         ylab("log2 p-value adjusted") +
         scale_colour_gradient(limits=c(0, 1), low="red", high="black") +
         theme_bw()
ggsave(file=paste0("volcanoplot", "~", variable_names[1], ".pdf") , width = 20, height = 20, units = "cm")

distsRL = dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hc <- hclust(distsRL)

pdf(paste0("Sample-to-sample distances~", variable_names[1], ".pdf") , height=10,width=10)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))
x = dev.off()











