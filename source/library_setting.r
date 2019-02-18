# .libPaths()



.libPaths( c( .libPaths(), "/BiO/apps/R-3.5.2/lib64/R/library/") )


myPath <- .libPaths()[2]


.libPaths(myPath)