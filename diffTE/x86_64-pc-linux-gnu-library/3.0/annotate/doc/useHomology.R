### R code from vignette source 'useHomology.Rnw'

###################################################
### code chunk number 1: loadAndListPackages
###################################################
library("hom.Hs.inp.db")
ls("package:hom.Hs.inp.db")


###################################################
### code chunk number 2: listMapping
###################################################
as.list(hom.Hs.inpMUSMU)[1:4]


###################################################
### code chunk number 3: mapGeneExample
###################################################
# load the organism annotation data for human
library(org.Hs.eg.db)

# get the entrex gene ID for gene symbol "MSX2"
mget("MSX2", org.Hs.egSYMBOL2EG)

# get the ensembl protein ID for the entrez gene ID "4488"
mget("4488", org.Hs.egENSEMBLPROT)

# use the inparanoid package to get the mouse gene that is considered 
# equivalent to ensembl protein ID "ENSP00000239243"
mget("ENSP00000239243", hom.Hs.inpMUSMU)

# load the organism annotation data for mouse
library(org.Mm.eg.db)

# get the entrez gene ID for Jackson labs ID "MGI:97169"
mget("MGI:97169", org.Mm.egMGI2EG)

# finally get the gene symbol for entrez gene ID "17702" from mouse
mget("17702", org.Mm.egSYMBOL)


###################################################
### code chunk number 4: seedPairExample
###################################################
mget("ENSP00000301011", hom.Hs.inpMUSMU)


###################################################
### code chunk number 5: seedPairExample2
###################################################
# make a connection to the human database
mycon <- hom.Hs.inp_dbconn()
# make a list of all the tables that are available in the DB
dbListTables(mycon)
# make a list of the columns in the table of interest
dbListFields(mycon, "mus_musculus")


###################################################
### code chunk number 6: seedPairExample3
###################################################
#make a query that will let us see which clust_id we need
sql <- "SELECT * FROM mus_musculus WHERE inp_id = 'ENSP00000301011';"
#retrieve the data
dataOut <- dbGetQuery(mycon, sql)
dataOut


###################################################
### code chunk number 7: seedPairExample4
###################################################
#make a query that will let us see all the data that is affiliated with a clust id
sql <- "SELECT * FROM mus_musculus WHERE clust_id = '1731';"
#retrieve the data
dataOut <- dbGetQuery(mycon, sql)
dataOut


###################################################
### code chunk number 8: useHomology.Rnw:257-258
###################################################
toLatex(sessionInfo())


