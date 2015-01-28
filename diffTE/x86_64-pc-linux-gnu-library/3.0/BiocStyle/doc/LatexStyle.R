### R code from vignette source 'LatexStyle.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: figureexample
###################################################
v = seq(0, 60i, length=1000)
plot(abs(v)*exp(v), type="l", col="Royalblue")


###################################################
### code chunk number 3: sessionInfo
###################################################
toLatex(sessionInfo())


