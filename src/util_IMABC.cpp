#include "util_IMABC.h"

// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Struct declarations
/*
struct Point{
    Point(std::vector<double> params, double pval, double dist, std::vector<double> pvals) :
    params(params), pval(pval), dist2_to_target(dist), pvals(pvals) {};

    std::vector<double> params;
    double pval;
    double dist2_to_target;
    std::vector<double> pvals;
};
*/
// Helper functions
double squared_dist_vec(std::vector<double> v1, std::vector<double> v2){
    if(v1.size() != v2.size()) {
        throw std::invalid_argument("Two vectors are not of the same length!");
    }

    size_t vsize = v1.size();
    double temp = 0.0;
    double squared_sum = 0.0;

    for(unsigned i=0;i<vsize; i++){
        squared_sum += (v1[i]-v2[i]) * (v1[i]-v2[i]);
    }

    return squared_sum;
}

/*const double log2pi = std::log(2.0 * M_PI);*/
arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd) {
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

    for (unsigned i=0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
        out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
    }

    if (logd == false) {
        out = exp(out);
    }
    return(out);
}

arma::vec arma_sum_dnorm(arma::mat &x,
                         std::vector<std::vector<double>>& mean_vec,
                         std::vector<NumericMatrix>& sigma_vec){

    if(mean_vec.size() != sigma_vec.size()){
        throw std::invalid_argument("Mean vector and sigma vector should have the same size!");
    }

    auto mean_it = mean_vec.begin();
    auto sigma_it = sigma_vec.begin();
    arma::vec sum_dnorm(x.n_rows);

    for(; mean_it != mean_vec.end(); mean_it++, sigma_it++){
        arma::rowvec m = as<arma::rowvec>(wrap(*mean_it));
        arma::mat sig = as<arma::mat>(*sigma_it);
        sum_dnorm += dmvnrm_arma(x, m, sig, false);
    }

    return sum_dnorm;
}

/*struct better_fit_target
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
};*/

std::vector<double> cal_pval_vec(std::vector<double> &obs,
                                 std::vector<double> &target, std::vector<double> &target_sd){
    std::vector<double> pvals(obs.size(), 0);

    for(unsigned i=0; i<obs.size(); i++){
        pvals[i] = 2*R::pnorm(obs[i], target[i], target_sd[i], (obs[i]<target[i]), false);
    }

    return pvals;
}

bool is_acceptable(std::vector<double> &pval_vec, std::vector<double> &alpha_vec){
    for(unsigned i=0; i<pval_vec.size(); i++){
        if(pval_vec[i] < alpha_vec[i]){
            return false;
        }
    }
    return true;
}

/*struct bad_point{
    bad_point(std::vector<double> alpha_vec) : alpha_vec(alpha_vec) {};

    bool operator()(Point p){
        return !is_acceptable(p.pvals, alpha_vec);
    }

    std::vector<double> alpha_vec;
};*/

// Print functions
void print_point_vec(const std::vector<Point> &vpoints){
    for(auto it = vpoints.begin(); it != vpoints.end(); it++){
        for(unsigned int i=0; i<(*it).params.size(); i++){
            Rcout << (*it).params[i] << " ";
        }
        Rcout << "\t";
        Rcout << "pval = " << (*it).pval << "\t";
        Rcout << "dist2_to_target = " << (*it).dist2_to_target << "\t";
        Rcout << "pvals = ";
        for(int i=0; i<(*it).pvals.size(); i++){
            Rcout << (*it).pvals[i] << " ";
        }
        Rcout << "\n";
    }
}
// End of Print functions

// Write to file function
void write_csv_vpoints(const std::vector<Point> &vpoints, std::string file_name,
                       int nparams, int ntargets){
    std::ofstream myfile;
    myfile.open(file_name);

    for(unsigned i=0; i<nparams; i++){
        myfile << "v" << i << ",";
    }
    myfile << "minp,";
    myfile << "dist,";
    for(unsigned i=0; i<ntargets; i++){
        myfile << "p" << i << ",";
    }
    myfile << "\n";

    for(auto it = vpoints.begin(); it != vpoints.end(); it++){
        for(unsigned i=0; i<(*it).params.size(); i++){
            myfile << (*it).params[i] << ",";
        }
        myfile << (*it).pval << ",";
        myfile << (*it).dist2_to_target << ",";
        for(unsigned i=0; i<(*it).pvals.size(); i++){
            myfile << (*it).pvals[i] << ",";
        }
        myfile << "\n";
    }
    myfile.close();
}
// End of Write to file function
// End of Helper functions

