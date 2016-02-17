#include <Rcpp.h>
#include <internal.h>
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
extern "C" const char* package_name = "multiStateFlowReliability";
namespace multistateTurnip
{
	SEXP crudeMC_graphNEL(SEXP graph, SEXP capacity, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices);
	SEXP crudeMC_igraph(SEXP graph, SEXP capacity, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices);
	SEXP crudeMC_graphAM(SEXP graph, SEXP capacity, SEXP n, SEXP threshold, SEXP seed, SEXP interestVertices);
}
R_CallMethodDef callMethods[] = 
{
	{"crudeMC_igraph", (DL_FUNC)&multistateTurnip::crudeMC_igraph, 6},
	{"crudeMC_graphNEL", (DL_FUNC)&multistateTurnip::crudeMC_graphNEL, 6},
	{"crudeMC_graphAM", (DL_FUNC)&multistateTurnip::crudeMC_graphAM, 6},
	{NULL, NULL, 0}
};
RcppExport void R_init_multiStateFlowReliability(DllInfo *info)
{
	std::vector<R_CallMethodDef> callMethodsVector;
	R_CallMethodDef* packageCallMethods = callMethods;
	while(packageCallMethods->name != NULL) packageCallMethods++;
	callMethodsVector.insert(callMethodsVector.begin(), callMethods, packageCallMethods);

	R_CallMethodDef* RcppStartCallMethods = Rcpp_get_call();
	R_CallMethodDef* RcppEndCallMethods = RcppStartCallMethods;
	while(RcppEndCallMethods->name != NULL) RcppEndCallMethods++;
	callMethodsVector.insert(callMethodsVector.end(), RcppStartCallMethods, RcppEndCallMethods);
	R_CallMethodDef blank = {NULL, NULL, 0};
	callMethodsVector.push_back(blank);

	R_registerRoutines(info, NULL, &(callMethodsVector[0]), NULL, NULL);
	init_Rcpp_cache();
}
