.onLoad <- function(libname, pkgname)
{
	library.dynam(package="multiStateFlowReliability", chname="multiStateFlowReliability", lib.loc = .libPaths())
}
