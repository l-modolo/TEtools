# TEtools

TEtools is composed of three tools (countTE, diffTE and PingPong) and their [galaxy website](http://getgalaxy.org/) wrapper.

## Install

To install TEtools in your galaxy instance
simply go to the `galaxy-dist/tools` folder of your galaxy instance folder and execute the following command:

```sh
git clone https://github.com/l-modolo/TEtools
```

Then you have to edit the file `galaxy-dist/config/tool_conf.xml` to add the TEtools to your tools menu by adding the following lines between the `<toolbox></toolbox>` tags.
```xml
<section id="TEtools" name="TEtools">
    <tool file="TEtools/countTE.xml" />
    <tool file="TEtools/diffTE.xml" />
    <!-- <tool file="TEtools/PingPong.xml" /> -->
</section>
```
At least you need to make sure that all the R dependecy for the `diffTE` tool are installed on yout system:

for `.deb` systems:
```sh
sudo apt-get install libxml2-dev
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

## countTE

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
countTE will group together the count of reads mapping on TE copies according to the column number defined in the second field.

### NGS data file

The NGS data set can be of two types: fastq sequence files or sam alignement files

#### fastq files
You can add any number of **fastq files** to be mapped on the fasta file, for paired-end data you must add the same number of paired fastq files.

When fastq files are provided, countTE computes an index of the fasta file and then maps the reads using [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) or [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
For smallRNA sequencing data we recommend to use [bowtie](http://bowtie-bio.sourceforge.net/index.shtml), which seems to perform better than [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
When using RNA sequencing data we recommend to use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and to specify the correct insert size used to build the library.

#### sam file
counTE outputs the sam alignment files corresponding to each fastq file (or pair of fastq files in the case of paired-end data) in the same order than these fastq files.
You can also directly use sam alignement files instead of fastq files to skip the mapping step of countTE.
This is useful when you want to compute a count table according to another column in the rosette file for example.

### output file
countTE reports a space delimited tabular text file of the read counts.

    - The first column corresponds to the rosette file column on which the read count was performed.
    - If more than one variable column was provided in the rosette file, they will be put after the first column.
    - The next following column(s), but the last, correspond to the number of mapping reads for each sample (fastq/sam files).
    - The last column corresponds to the total of these counts.

### counTE.ini file
Some options for countTE are not available through the command line options. They are defined in a counTE.ini file. This file contains options like the size of a siRNA, the number of threads to use or the path of the different programs called by counTE. By defaut, this countTE.ini file is created with default options, if not found in the same directory as the file countTE.py.

## diffTE

diffTE performs a differential expression analysis on the counTE output file using [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html).

This tool produces an HTML output with clickable images, allowing to download PDF files and link to a table of differentially expressed TEs.


## PingPong
in development

