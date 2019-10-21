#ifndef UTIL_IMABC_H
#define UTIL_IMABC_H

// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Struct declarations
struct Point{
    Point(std::vector<double> params, double pval, double dist, std::vector<double> pvals) :
    params(params), pval(pval), dist2_to_target(dist), pvals(pvals) {};

    std::vector<double> params;
    double pval;
    double dist2_to_target;
    std::vector<double> pvals;
};

// Helper functions
double squared_dist_vec(std::vector<double> v1, std::vector<double> v2);

const double log2pi = std::log(2.0 * M_PI);

arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma);

arma::vec arma_sum_dnorm(arma::mat &x,
                         std::vector<std::vector<double>>& mean_vec,
                         std::vector<NumericMatrix>& sigma_vec);

struct better_fit_target
{
    inline bool operator() (const Point& p1, const Point& p2){
        return (p1.pval == p2.pval ? p1.dist2_to_target < p2.dist2_to_target : p1.pval > p2.pval);
    }
};

struct less_dist_to_ref
{
    less_dist_to_ref(Point ref) : ref(ref) {};
    inline bool operator() (const Point& p1, const Point& p2){
        return (squared_dist_vec(ref.params, p1.params) < squared_dist_vec(ref.params, p2.params));
    }
    Point ref;
};

std::vector<double> cal_pval_vec(std::vector<double> &obs,
                                 std::vector<double> &target, std::vector<double> &target_sd);

bool is_acceptable(std::vector<double> &pval_vec, std::vector<double> &alpha_vec);

struct bad_point{
    bad_point(std::vector<double> alpha_vec) : alpha_vec(alpha_vec) {};

    bool operator()(Point p){
        return !is_acceptable(p.pvals, alpha_vec);
    }

    std::vector<double> alpha_vec;
};

// Print functions
void print_point_vec(const std::vector<Point> &vpoints);
// End of Print functions

// Write to file function
void write_csv_vpoints(const std::vector<Point> &vpoints, std::string file_name,
                       int nparams, int ntargets);
// End of Write to file function
// End of Helper functions

#endif
