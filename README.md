# TEtools

TEtools is composed of three tools (TEcount, TEdiff and PingPong) and their [galaxy website](http://getgalaxy.org/) wrapper.

## Install

To install TEtools in your galaxy instance
simply go to the `galaxy-dist/tools` folder of your galaxy instance folder and execute the following command:

```sh
git clone https://github.com/l-modolo/TEtools
```

Then you have to edit the file `galaxy-dist/config/tool_conf.xml` to add the TEtools to your tools menu by adding the following lines between the `<toolbox></toolbox>` tags.
```xml
<section id="TEtools" name="TEtools">
    <tool file="TEtools/TEcount.xml" />
    <tool file="TEtools/TEdiff.xml" />
    <!-- <tool file="TEtools/PingPong.xml" /> -->
</section>
```
At least you need to make sure that all the R dependecy for the `TEdiff` tool are installed on yout system:

for `.deb` systems:
```sh
sudo apt-get install libxml2-dev python-configparser
```
for `.rpm`systems:
```sh
yum install libxml2-devel glibc
```

then execute the following commands in an [R](http://cran.r-project.org/) console:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2", dep=T)
biocLite("gplots", dep=T)
biocLite("ggplot2", dep=T)
biocLite("RColorBrewer", dep=T)
```

Restart your galaxy server and enjoy TEtools !

## TEcount

This program computes a count table file using NGS data file(s), a fasta file containing TE copy sequences and a rosette file.

### rosette file
The rosette file contains at least 2 columns. The first column corresponds to the names of the TE copies as in the fasta file, and the second column corresponds to a variable associated to these TE copy names on which we want to compute the counts (for example TE familly).

For example, we can write the following rosette file:
```
2L|(3071416..3071503,3071708..3071841)|DNA/P|PROTOP   PROTOP
2L|(5363113..5363154,5363819..5363952)|DNA/P|PROTOP   PROTOP
2L|c(9889960..9890093,9890313..9890400)|DNA/P|PROTOP  PROTOP
2L|(20948958..20949699)|DNA/RC|DNAREP1_DM             DNAREP1_DM
2L|c(20958914..20959207)|DNA/RC|DNAREP1_DM            DNAREP1_DM
2L|c(20966385..20966456)|DNA/RC|DNAREP1_DM            DNAREP1_DM
2L|(20976274..20976387)|DNA/RC|DNAREP1_DM             DNAREP1_DM
```

This will allow to count reads mapping on the `PROTOP` and the `DNAREP1_DM` elements.
The rosette file can contain more TE copy names than there is sequences in the fasta file, but we cannot map a read on a TE copy not present in the fasta file.
The fasta file can contain copies not present in the rosette file, but reads mapping on these copies will be ignored.

The rosette file can contain as many variable columns as necessary.
TEcount will group together the count of reads mapping on TE copies according to the column number defined in the second field.

### NGS data file

The NGS data set can be of two types: fastq sequence files or sam alignement files

#### fastq files
You can add any number of **fastq files** to be mapped on the fasta file, for paired-end data you must add the same number of paired fastq files.

When fastq files are provided, TEcount computes an index of the fasta file and then maps the reads using [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) or [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
For smallRNA sequencing data we recommend to use [bowtie](http://bowtie-bio.sourceforge.net/index.shtml), which seems to perform better than [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
When using RNA sequencing data we recommend to use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and to specify the correct insert size used to build the library.

#### sam file
TEcount outputs the sam alignment files corresponding to each fastq file (or pair of fastq files in the case of paired-end data) in the same order than these fastq files.
You can also directly use sam alignement files instead of fastq files to skip the mapping step of TEcount.
This is useful when you want to compute a count table according to another column in the rosette file for example.

### output file
TEcount reports a space delimited tabular text file of the read counts.

    - The first column corresponds to the rosette file column on which the read count was performed.
    - If more than one variable column was provided in the rosette file, they will be put after the first column.
    - The next following column(s), but the last, correspond to the number of mapping reads for each sample (fastq/sam files).
    - The last column corresponds to the total of these counts.

### TEcount.ini file
Some options for TEcount are not available through the command line options. They are defined in a TEcount.ini file. This file contains options like the size of a siRNA, the number of threads to use or the path of the different programs called by TEcount. By defaut, this TEcount.ini file is created with default options, if not found in the same directory as the file TEcount.py.

## TEdiff

TEdiff performs a differential expression analysis on the TEcount output file using [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html).

This tool produces an HTML output with clickable images, allowing to download PDF files and link to a table of differentially expressed TEs.


## PingPong
in development


# TEtools command line usage
As TEtools was developed for a galaxy interface, the different parts of the pipeline can be manually called from a command line interface.
The installation procedure without a galaxy instance for an exclusive command line usage is exactly the same, you just can clone the TEtools repository in a place of your choosing instead of the `galaxy-dist/tools` folder.

## TEcount
```sh
TEcount.py -rosette [$]ROSETTE_FILE] -column [$]COUNT_COLUMN] -TE_fasta [FASTA_FILE] -count [OUTPUT_FILE] -RNA [FASTQ_FILE1 FASTQ_FILE2 ... FASTQ_FILEN]
```
For RNASeq data you can add the option `-bowtie2` and to run the UrQt quality trimming software before the mapping, you can add the option `-QC`.
For paired-end data the second list of fastq files must be entered after the option `-RNApair [FASTQ_FILE1 FASTQ_FILE2 ... FASTQ_FILEN]` and the insert size must be specified with the option `-insert [SIZE]`.
Like with galaxy you can directly use sam files instead of fastq files by replacing the option `-RNA` by `-sam [SAM_FILE1 SAM_FILE2 ... SAM_FILEN]`.
An alternative count file for reads of size 21pb can be computed with the option `-siRNA [OUTPUT_SIRNA_FILE]`

## TEdiff
```sh
Rscript TEdiff.R --args --FDR_level=[FDR_LEVEL] --count_column=[COUNT_COLUMN] --count_file=\"[COUNT_FILE]\" experiment_formula=\"[EXPERIMENT_FORMULA]\" --sample_names=\"[SAMPLE_NAMES]\" --outdir=\"[OUTPUT_HTML_FOLDER]\" --htmlfile=\"[OUTPUT_HTML_FILE]\"
```

With:
+ `[FDR_LEVEL]` the FDR threshold for the analysis
+ `[COUNT_FILE]` the count file computed with TEcount
+ `[COUNT_COLUMN]` the number of the first column with read counts (ex: 2 if the rosette file as only one variable)
+ `[EXPERIMENT_FORMULA]` the formula corresponding to the names (example: `sample_name:replicat:condition)
+ `[SAMPLE_NAMES]` a comma separated list of names of the different samples in the same order as in the count file (example for n samples with 2 replicats: `sample_name_1:replicat_1,sample_name_1:replicat_2,...,sample_name_n:replicat_1,sample_name_n:replicat_2)`


