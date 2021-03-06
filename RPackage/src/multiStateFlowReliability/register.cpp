#include <Rcpp.h>
#include <internal.h>
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
extern "C" const char* package_name = "multiStateFlowReliability";
#include "RPackage/crudeMC.h"
#include "RPackage/PMC.h"
#include "RPackage/turnip.h"
namespace multistateTurnip
{
	SEXP generalisedSplittingFixedEffort_igraph(SEXP graph, SEXP distributions, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp);
	SEXP generalisedSplittingFixedEffort_graphNEL(SEXP graph, SEXP distributions, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp);
	SEXP generalisedSplittingFixedEffort_graphAM(SEXP graph, SEXP distributions, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp);
	SEXP generalisedSplittingFixedFactors_igraph(SEXP graph, SEXP distributions, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp);
	SEXP generalisedSplittingFixedFactors_graphNEL(SEXP graph, SEXP distributions, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp);
	SEXP generalisedSplittingFixedFactors_graphAM(SEXP graph, SEXP distributions, SEXP n, SEXP levels_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp);
	SEXP generalisedSplittingFixedFactorsEvolution_igraph(SEXP graph, SEXP distributions, SEXP n, SEXP times_sexp, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp);
	SEXP generalisedSplittingFixedFactorsEvolution_graphNEL(SEXP graph, SEXP distributions, SEXP n, SEXP times_sexp, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp);
	SEXP generalisedSplittingFixedFactorsEvolution_graphAM(SEXP graph, SEXP distributions, SEXP n, SEXP times_sexp, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP factors_sexp);
	SEXP generalisedSplittingAdaptiveEvolution_igraph(SEXP graph, SEXP distributions, SEXP n, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP fraction_sexp);
	SEXP generalisedSplittingAdaptiveEvolution_graphNEL(SEXP graph, SEXP distributions, SEXP n, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP fraction_sexp);
	SEXP generalisedSplittingAdaptiveEvolution_graphAM(SEXP graph, SEXP distributions, SEXP n, SEXP level_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, SEXP verbose_sexp, SEXP fraction_sexp);

}
R_CallMethodDef callMethods[] = 
{
	{"crudeMC_igraph", (DL_FUNC)&multistateTurnip::crudeMC_igraph, 7},
	{"crudeMC_graphNEL", (DL_FUNC)&multistateTurnip::crudeMC_graphNEL, 7},
	{"crudeMC_graphAM", (DL_FUNC)&multistateTurnip::crudeMC_graphAM, 7},
	{"pmc_igraph", (DL_FUNC)&multistateTurnip::pmc_igraph, 8},
	{"pmc_graphNEL", (DL_FUNC)&multistateTurnip::pmc_graphNEL, 8},
	{"pmc_graphAM", (DL_FUNC)&multistateTurnip::pmc_graphAM, 8},
	{"turnip_igraph", (DL_FUNC)&multistateTurnip::turnip_igraph, 9},
	{"turnip_graphNEL", (DL_FUNC)&multistateTurnip::turnip_graphNEL, 9},
	{"turnip_graphAM", (DL_FUNC)&multistateTurnip::turnip_graphAM, 9},
	{"generalisedSplittingFixedEffort_igraph", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedEffort_igraph, 7},
	{"generalisedSplittingFixedEffort_graphNEL", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedEffort_graphNEL, 7},
	{"generalisedSplittingFixedEffort_graphAM", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedEffort_graphAM, 7},
	{"generalisedSplittingFixedFactors_igraph", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedFactors_igraph, 8},
	{"generalisedSplittingFixedFactors_graphNEL", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedFactors_graphNEL, 8},
	{"generalisedSplittingFixedFactors_graphAM", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedFactors_graphAM, 8},
	{"generalisedSplittingFixedFactorsEvolution_igraph", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedFactorsEvolution_igraph, 9},
	{"generalisedSplittingFixedFactorsEvolution_graphNEL", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedFactorsEvolution_graphNEL, 9},
	{"generalisedSplittingFixedFactorsEvolution_graphAM", (DL_FUNC)&multistateTurnip::generalisedSplittingFixedFactorsEvolution_graphAM, 9},
	{"generalisedSplittingAdaptiveEvolution_igraph", (DL_FUNC)&multistateTurnip::generalisedSplittingAdaptiveEvolution_igraph, 8},
	{"generalisedSplittingAdaptiveEvolution_graphNEL", (DL_FUNC)&multistateTurnip::generalisedSplittingAdaptiveEvolution_graphNEL, 8},
	{"generalisedSplittingAdaptiveEvolution_graphAM", (DL_FUNC)&multistateTurnip::generalisedSplittingAdaptiveEvolution_graphAM, 8},
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
