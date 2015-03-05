TEtools
=======

TEtools is composed of three tools and their [galaxy website](http://getgalaxy.org/) wrapper.

Install
-------

To install TEtools in your galaxy instance
simply go to the `galaxy-dist/tools` folder of your galaxy instance folder and execute the following command:

```sh
git clone https://github.com/l-modolo/TEtools
```

Then you have to edit the file `galaxy-dist/config/tool_conf.xml` to add the TEtools to your tools menu by adding the following lines between the `<toolbox></toolbox>` tags.
```sh
<section id="TEtools" name="TEtools">
    <tool file="TEtools/countTE.xml" />
    <tool file="TEtools/diffTE.xml" />
    <tool file="TEtools/PingPong.xml" />
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

then execute the following commands in an `R` console:
```R
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2", dep=T)
biocLite("gplots", dep=T)
biocLite("ggplot2", dep=T)
biocLite("RColorBrewer", dep=T)
```

Restart your galaxy server and enjoy TEtools !

