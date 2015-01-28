## ----------------------------------------------------------------------
## This file contains the steps that were used 
## 1.+2. to prepare the example datasets "lymphoma" and "kidney"
##       of the package in 2002
## 3.+4. update them from exprSet to ExpressionSet on 9 Feb 2007.
## ----------------------------------------------------------------------
dataout   = "../../data"
lym = file.path(dataout, 'lymphoma.RData')
kid = file.path(dataout, 'kidney.RData')

library("Biobase")

## ------------------------------------------------------------
## 1. create lymphoma
## ------------------------------------------------------------
samples   = read.delim("lymphomasamples.txt", as.is=T)
datain    = "/home/whuber/h/VSN/alizadeh"

## CH1  = Cy3 = green =  reference
## CH2  = Cy5 = red   =  sample of interest

nrspots   = 9216
nrsamples = nrow(samples)
qua       = matrix(NA, nrow=nrspots, ncol=2*nrsamples)
pd        = data.frame(name    = I(character(2*nrsamples)),
                       sample  = I(character(2*nrsamples)))
  
for (i in 1:nrsamples) {
  filename = paste(samples$name[i], 'rex.DAT', sep='')
  dat = read.delim(file.path(datain, filename))
  qua[,2*i-1] = dat$CH1I - dat$CH1B
  qua[,2*i]   = dat$CH2I - dat$CH2B
  pd$name[(2*i-1):(2*i)]    = samples$name[i]
  pd$sample[2*i-1]  = "reference"
  pd$sample[2*i]    = samples$sampleid[i]
}
colnames(qua) = pd$sample

lymphomaPhenoData <- new("AnnotatedDataFrame")
pData(lymphomaPhenoData) <- pd
varLabels(lymphomaPhenoData) <- list(name="Name of the Chip", sample="Sample")

lymphoma = new("ExpressionSet",
    exprs = qua,
    phenoData = lymphomaPhenoData)
    
save(lymphoma, file=lym, compress=TRUE)

## ------------------------------------------------------------
## 2. create kidney
## ------------------------------------------------------------
datain = "/net/herkules/raid4/home/whuber/Kidney2"
thehyb = 90
load(file.path(datain, "squa.Rdata"))

dat = (squa[, c("fg.green", "fg.red"), thehyb]
      -squa[, c("bg.green", "bg.red"), thehyb])
rownames(dat) = NULL
colnames(dat) = c("green", "red")

kidneyPhenoData <- new("AnnotatedDataFrame")
pData(kidneyPhenoData) <- data.frame(channel = c("green", "red"))
varLabels(kidneyPhenoData) <- list(channel="green: 532 nm, dye=Cy3; red: 635 nm, dye=Cy5")

kidney = new("ExpressionSet",
  exprs = dat,
  phenoData = kidneyPhenoData)

save(kidney, file=kid, compress=TRUE)




## ------------------------------------------------------------
## 3. create lymphoma
## ------------------------------------------------------------
load(lym)
pd = pData(lymphoma)
ex = exprs(lymphoma)

pd$dye = factor(paste("Cy", rep(c(3L, 5L), 8), sep=""))

rownames(pd) = colnames(ex) = paste(pd$name, pd$sample, sep=".")

vmd = data.frame(labelDescription=I(c("Array ID", "Sample", "Dye")))
rownames(vmd) = colnames(pd)
ad = new("AnnotatedDataFrame", data=pd, varMetadata=vmd)

lymphoma = new("ExpressionSet", exprs=ex, phenoData=ad)

save(lymphoma, file=lym)


## ------------------------------------------------------------
## 4. update kidney
## ------------------------------------------------------------
load(kid)

pd =  pData(kidney)
ex =  exprs(kidney)
rownames(pd) = colnames(ex)

vmd = data.frame(labelDescription=I("The scanner channel Cy3 or Cy5"))
rownames(vmd) = colnames(pd)
ad = new("AnnotatedDataFrame", data=pd, varMetadata=vmd)
    
kidney = new("ExpressionSet", exprs=ex, phenoData=ad)

save(kidney, file=kid)
