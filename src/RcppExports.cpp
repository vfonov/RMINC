// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// anat_summary
List anat_summary(CharacterVector filenames, IntegerVector atlas, std::string method);
RcppExport SEXP RMINC_anat_summary(SEXP filenamesSEXP, SEXP atlasSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filenames(filenamesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type atlas(atlasSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(anat_summary(filenames, atlas, method));
    return rcpp_result_gen;
END_RCPP
}
// count_labels
NumericMatrix count_labels(CharacterVector filenames);
RcppExport SEXP RMINC_count_labels(SEXP filenamesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filenames(filenamesSEXP);
    rcpp_result_gen = Rcpp::wrap(count_labels(filenames));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_minc_apply
List rcpp_minc_apply(CharacterVector filenames, bool use_mask, CharacterVector mask, double mask_lower_val, double mask_upper_val, RObject value_for_mask, bool filter_masked, NumericVector slab_sizes, Function fun, List args);
RcppExport SEXP RMINC_rcpp_minc_apply(SEXP filenamesSEXP, SEXP use_maskSEXP, SEXP maskSEXP, SEXP mask_lower_valSEXP, SEXP mask_upper_valSEXP, SEXP value_for_maskSEXP, SEXP filter_maskedSEXP, SEXP slab_sizesSEXP, SEXP funSEXP, SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type filenames(filenamesSEXP);
    Rcpp::traits::input_parameter< bool >::type use_mask(use_maskSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< double >::type mask_lower_val(mask_lower_valSEXP);
    Rcpp::traits::input_parameter< double >::type mask_upper_val(mask_upper_valSEXP);
    Rcpp::traits::input_parameter< RObject >::type value_for_mask(value_for_maskSEXP);
    Rcpp::traits::input_parameter< bool >::type filter_masked(filter_maskedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type slab_sizes(slab_sizesSEXP);
    Rcpp::traits::input_parameter< Function >::type fun(funSEXP);
    Rcpp::traits::input_parameter< List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_minc_apply(filenames, use_mask, mask, mask_lower_val, mask_upper_val, value_for_mask, filter_masked, slab_sizes, fun, args));
    return rcpp_result_gen;
END_RCPP
}
// graph_tfce_wqu
std::vector<double> graph_tfce_wqu(std::vector<double> map, std::vector<std::vector<int> > adjacencies, double E, double H, int nsteps, std::vector<double> weights);
RcppExport SEXP RMINC_graph_tfce_wqu(SEXP mapSEXP, SEXP adjacenciesSEXP, SEXP ESEXP, SEXP HSEXP, SEXP nstepsSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type map(mapSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type adjacencies(adjacenciesSEXP);
    Rcpp::traits::input_parameter< double >::type E(ESEXP);
    Rcpp::traits::input_parameter< double >::type H(HSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_tfce_wqu(map, adjacencies, E, H, nsteps, weights));
    return rcpp_result_gen;
END_RCPP
}
// graph_tfce
std::vector<double> graph_tfce(std::vector<double> map, std::vector<std::vector<int> > adjacencies, double E, double H, int nsteps, std::vector<double> weights);
RcppExport SEXP RMINC_graph_tfce(SEXP mapSEXP, SEXP adjacenciesSEXP, SEXP ESEXP, SEXP HSEXP, SEXP nstepsSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type map(mapSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type adjacencies(adjacenciesSEXP);
    Rcpp::traits::input_parameter< double >::type E(ESEXP);
    Rcpp::traits::input_parameter< double >::type H(HSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_tfce(map, adjacencies, E, H, nsteps, weights));
    return rcpp_result_gen;
END_RCPP
}
// coords2ind
int coords2ind(int i, int j, int k, int d1, int d2, int d3);
RcppExport SEXP RMINC_coords2ind(SEXP iSEXP, SEXP jSEXP, SEXP kSEXP, SEXP d1SEXP, SEXP d2SEXP, SEXP d3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< int >::type d2(d2SEXP);
    Rcpp::traits::input_parameter< int >::type d3(d3SEXP);
    rcpp_result_gen = Rcpp::wrap(coords2ind(i, j, k, d1, d2, d3));
    return rcpp_result_gen;
END_RCPP
}
// ind2coords
std::vector<int> ind2coords(int v, int d1, int d2, int d3);
RcppExport SEXP RMINC_ind2coords(SEXP vSEXP, SEXP d1SEXP, SEXP d2SEXP, SEXP d3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< int >::type d2(d2SEXP);
    Rcpp::traits::input_parameter< int >::type d3(d3SEXP);
    rcpp_result_gen = Rcpp::wrap(ind2coords(v, d1, d2, d3));
    return rcpp_result_gen;
END_RCPP
}
// neighbour_list
std::vector<std::vector<int> > neighbour_list(double x, double y, double z, int n);
RcppExport SEXP RMINC_neighbour_list(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(neighbour_list(x, y, z, n));
    return rcpp_result_gen;
END_RCPP
}
// mesh_area
std::vector<double> mesh_area(std::vector<double> vertices, std::vector<double> triangles);
RcppExport SEXP RMINC_mesh_area(SEXP verticesSEXP, SEXP trianglesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type vertices(verticesSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type triangles(trianglesSEXP);
    rcpp_result_gen = Rcpp::wrap(mesh_area(vertices, triangles));
    return rcpp_result_gen;
END_RCPP
}
