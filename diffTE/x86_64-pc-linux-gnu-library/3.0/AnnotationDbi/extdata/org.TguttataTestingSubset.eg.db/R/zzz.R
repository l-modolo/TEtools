datacache <- new.env(hash=TRUE, parent=emptyenv())

org.TguttataTestingSubset.eg <- function() showQCData("org.TguttataTestingSubset.eg", datacache)
org.TguttataTestingSubset.eg_dbconn <- function() dbconn(datacache)
org.TguttataTestingSubset.eg_dbfile <- function() dbfile(datacache)
org.TguttataTestingSubset.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.TguttataTestingSubset.eg_dbInfo <- function() dbInfo(datacache)

org.TguttataTestingSubset.egORGANISM <- "Taeniopygia guttataTestingSubset"

.onLoad <- function(libname, pkgname)
{
    require("methods", quietly=TRUE)
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.TguttataTestingSubset.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.TguttataTestingSubset.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.TguttataTestingSubset.eg_dbconn())
}

