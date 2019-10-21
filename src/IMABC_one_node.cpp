// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "util_IMABC.h"
#include "sim_one_node.h"
using namespace Rcpp;

class IMABC_one_node {
public:
    IMABC_one_node(size_t size_theta,
          std::vector<double> target,
          std::vector<double> target_sd,
          std::vector<double> alpha_start,
          std::vector<double> alpha_goal,
          int N0, int Nc, int Ngoal, int B, int LIM1, int LIM2, int LIM3) :
        size_theta(size_theta),
        target(target), target_sd(target_sd),
        alpha_vec(alpha_start), alpha_goal(alpha_goal),
        N0(N0), Nc(Nc), Ngoal(Ngoal), B(B), LIM1(LIM1), LIM2(LIM2), LIM3(LIM3) {
        ret = std::vector<Point>();
        theta_vec = std::vector<std::vector<double>>();
        sigma_vec = std::vector<NumericMatrix>();
        running_Nt = 0;
        n_selected_pts = 0;
        ESS = 0.0;
    };

    double d_prior(std::vector<double> x){
        double d_prior_val = R::punif(x[0], 60, 85, false, false)
                        + R::punif(x[1], 0, 5, false, false)
                        + R::punif(x[2], 5, 10, false, false)
                        + R::punif(x[3], 40, 85, false, false)
                        + R::punif(x[4], 5, 30, false, false)
                        + R::punif(x[5], 0, 10, false, false)
                        + R::punif(x[6], 0, 10, false, false);

        return d_prior_val;
    }

    // Microsimulation function
    std::vector<double> fsim(std::vector<double> x){
        Sim_one_node s(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
        return s.read_csv_and_schedule_cancer();
    }

    double var_unif(double a, double b){
        return (b-a)*(b-a)/12;
    }

    void step1(){
        std::vector<Point>().swap(ret);
        std::vector<double> params(size_theta);
        double min_p, dist2_to_target;

        for(unsigned it=0; it<N0; it++){
            /*
            for(unsigned i=0; i<size_theta; ++i){
                params[i] = R::rnorm(3, 1);
            }
            */
            params[0] = R::runif(60, 85);
            params[1] = R::runif(0, 5);
            params[2] = R::runif(5, 10);
            params[3] = R::runif(40, 85);
            params[4] = R::runif(5, 30);
            params[5] = R::runif(0, 10);
            params[6] = R::runif(0, 10);


            std::vector<double> curr_sim_res = fsim(params);
            std::vector<double> pval_vec = cal_pval_vec(curr_sim_res, target, target_sd);
            min_p = *(std::min_element(pval_vec.begin(), pval_vec.end()));
            dist2_to_target = squared_dist_vec(curr_sim_res, target);
            ret.push_back(Point(params, min_p, dist2_to_target, pval_vec));
        }

        running_Nt += N0;
    }

    void step2(){
        std::sort(ret.begin(), ret.end(), better_fit_target());
        Rcout << "Done step 2\n";
    }

    void step3(){
        std::vector<Point> ret_copy(ret);
        std::vector<Point> new_points = std::vector<Point>();
        std::vector<double> curr_sim_res;
        NumericMatrix sig_mat;
        NumericMatrix new_vecs;
        unsigned i;
        std::vector<Point>::iterator it, it_copy, it_begin, it_Nc;

        Environment mvrnorm_env("package:MASS");
        Function mvrnorm = mvrnorm_env["mvrnorm"];
        Function var("var");
        n_selected_pts = ret.size();

        if(!(n_selected_pts > LIM1*size_theta)){
            sig_mat = NumericMatrix::diag(size_theta, 1);
            sig_mat(0, 0) = var_unif(60, 85)/4;
            sig_mat(1, 1) = var_unif(0, 5)/4;
            sig_mat(2, 2) = var_unif(5, 10)/4;
            sig_mat(3, 3) = var_unif(40, 85)/4;
            sig_mat(4, 4) = var_unif(5, 30)/4;
            sig_mat(5, 5) = var_unif(0, 10)/4;
            sig_mat(6, 6) = var_unif(0, 10)/4;

        } else if((n_selected_pts > LIM1*size_theta) && (!(n_selected_pts > LIM2*size_theta))) {
            NumericMatrix param_mat(ret.size(), size_theta);
            for(i=0, it=ret.begin(); it!=ret.end(); it++, i++){
                NumericVector curr_params = wrap((*it).params);
                param_mat.row(i) = curr_params;
            }
            sig_mat = var(param_mat);
        } else {
        }
        Rcout << "In step 3\n";
        std::vector<double> m, curr_new_vec, pval_vec;
        double min_pval, dist2_to_target;

        it_begin=ret.begin();
        it_Nc= ret.size()<Nc ? ret.end() : ret.begin()+Nc;

        for(std::vector<Point>::iterator it = it_begin; it < it_Nc; it++){
            m = (*it).params;
            curr_sim_res = fsim(m);
            pval_vec = cal_pval_vec(curr_sim_res, target, target_sd);
            min_pval = *(std::min_element(pval_vec.begin(), pval_vec.end()));
            dist2_to_target = squared_dist_vec(curr_sim_res, target);
            if(is_acceptable(pval_vec, alpha_vec)){
                new_points.push_back(Point(m, min_pval, dist2_to_target, pval_vec));
            }

            if(n_selected_pts > LIM2*size_theta){
                auto ref_point = (*it);
                int n_nearest = LIM2 * size_theta;

                std::sort(ret_copy.begin(), ret_copy.end(), less_dist_to_ref(ref_point));

                NumericMatrix param_mat(n_nearest, size_theta);

                for(i=0, it_copy=ret_copy.begin(); it_copy!=ret_copy.begin() + n_nearest; it_copy++, i++){
                    NumericVector curr_params = wrap((*it_copy).params);
                    param_mat.row(i) = curr_params;
                }

                sig_mat = var(param_mat);
            }

            Rcout << "In step 3 2\n";

            if(B==1){
                new_vecs = NumericMatrix(1, size_theta);
                NumericVector temp = mvrnorm(B, m, sig_mat);
                new_vecs.row(0) = temp;
            } else {
                new_vecs = mvrnorm(B, m, sig_mat);
            }
            running_Nt += B;

            //Store the theta vec and sigma vec for the centre points
            theta_vec.push_back(m);
            sigma_vec.push_back(sig_mat);

            //Calculate the pval_vec and see if the point is acceptable
            for(i=0; i<B; i++){
                NumericVector temp_new_vec = new_vecs.row(i);
                curr_new_vec = as<std::vector<double>> (temp_new_vec);
                if(curr_new_vec[0]<60 || curr_new_vec[0]>85 ||
                   curr_new_vec[1]<0 || curr_new_vec[1]>5 ||
                   curr_new_vec[2]<5 || curr_new_vec[2]>10 ||
                   curr_new_vec[3]<40 || curr_new_vec[3]>85 ||
                   curr_new_vec[4]<5 || curr_new_vec[4]>30 ||
                   curr_new_vec[5]<0 || curr_new_vec[5]>10 ||
                   curr_new_vec[6]<0 || curr_new_vec[6]>10
                ) {
                    continue;
                }
                curr_sim_res = fsim(curr_new_vec);
                pval_vec = cal_pval_vec(curr_sim_res, target, target_sd);
                min_pval = *(std::min_element(pval_vec.begin(), pval_vec.end()));
                dist2_to_target = squared_dist_vec(curr_sim_res, target);

                if(is_acceptable(pval_vec, alpha_vec)){
                    new_points.push_back(Point(curr_new_vec, min_pval, dist2_to_target, pval_vec));
                }
            }
        }
        ret.erase(it_begin, it_Nc);
        ret.insert(ret.end(), new_points.begin(), new_points.end());
        n_selected_pts = ret.size();

        if(is_trial) Rcout << "running_Nt = " << running_Nt << "\n";
        if(is_trial) Rcout << "number of selected points: " << n_selected_pts << "\n";
    }

    // Calculate ESS
    void step4(){
        //double curr_d_prior, curr_dmixture;
        //double curr_weight = 1;
        arma::mat params_mat(ret.size(), size_theta);
        arma::vec dprior(ret.size());
        std::vector<Point>::iterator it;
        unsigned i;

        for(it = ret.begin(), i=0; it != ret.end(); it++, i++){
            std::vector<double> curr_params = (*it).params;
            dprior[i] = d_prior(curr_params);
            params_mat.row(i) = as<arma::vec>(wrap(curr_params)).t();
        }

        //Rf_PrintValue(wrap(dprior));

        arma::vec sum_dmixture = arma_sum_dnorm(params_mat, theta_vec, sigma_vec);
        arma::vec dmixture = (sum_dmixture*B + dprior * N0)/(B*theta_vec.size() + N0);
        arma::vec weight = dprior/dmixture;
        ESS = 1.0/arma::accu(arma::pow(weight, 2));
        Rcout << "ESS = " << ESS << "\n";
        write_csv_vpoints(ret, "myparamsonenode" + std::to_string(save_i++) + ".csv", size_theta, target.size());
    }

    void IMABC_main(){
        step1();
        write_csv_vpoints(ret, "myparamsonenode" + std::to_string(save_i++) + ".csv", size_theta, target.size());

        while(alpha_vec < alpha_goal){
            if(n_selected_pts > LIM3 * size_theta){
                std::sort(ret.begin(), ret.end(), better_fit_target());
                int median_idx = int(n_selected_pts/2);
                std::vector<Point>::iterator it = ret.begin() + median_idx;
                std::vector<double>pvals_median = (*it).pvals;
                for(unsigned i=0; i<alpha_goal.size(); i++){
                    alpha_vec[i] = std::min(pvals_median[i], alpha_goal[i]);
                }
                if(is_trial) Rcout << n_selected_pts << "\n";
                int n_preserved = int(n_selected_pts/4);
                ret.erase(std::remove_if(ret.begin() + n_preserved, ret.end(), bad_point(alpha_vec)),
                          ret.end());
                n_selected_pts = ret.size();
                if(is_trial) Rcout << n_selected_pts << "\n";
            } else {
                step2();
                step3();
                write_csv_vpoints(ret, "myparamsonenode" + std::to_string(save_i++) + ".csv", size_theta, target.size());
            }
        }

        step4();

        while(ESS < Ngoal){
            step2();
            step3();
            step4();
        }

    }

    // data
    size_t size_theta;
    std::vector<double> target;
    std::vector<double> target_sd;
    std::vector<double> alpha_vec;
    std::vector<double> alpha_goal;
    unsigned N0, Nc, Ngoal, B, LIM1, LIM2, LIM3, running_Nt=0, n_selected_pts=0;
    double ESS = 0.0;

    std::vector<Point> ret;
    std::vector<std::vector<double>> theta_vec;
    std::vector<NumericMatrix> sigma_vec;

    int save_i = 0;
    bool is_trial = true;
};

//' @export
//[[Rcpp::export]]
void create_IMABC_one_node(unsigned N0, unsigned Nc, unsigned Ngoal, unsigned B,
                  NumericVector target_sd_in,
                  NumericVector alpha_start_in,
                  NumericVector alpha_goal_in){
    size_t size_theta = 7;
    //int N0=1e5, Nc=20, Ngoal=2, B=2000,
    int LIM1=5, LIM2=25, LIM3=50;
    std::vector<double>target{0.0000000, 0.1826170, 0.4079510, 0.6388139, 1.0000000,
                              0.0000000, 0.2883196, 0.4851360, 0.5900274, 1.0000000,
                              0.0000000, 0.2769945, 0.5669079, 0.5669079, 1.0000000};
    /*std::vector<double>target_sd{0.1, 0.11, 0.12, 0.13, 0.14,
                                 0.1, 0.11, 0.12, 0.13, 0.14,
                                 1, 1, 1, 1, 1};

    std::vector<double>alpha_start(15, 0.0000001);
    std::vector<double>alpha_goal(15, 0.001);
    */
    std::vector<double>target_sd = as<std::vector<double>>(target_sd_in);
    std::vector<double>alpha_start = as<std::vector<double>>(alpha_start_in);
    std::vector<double>alpha_goal = as<std::vector<double>>(alpha_goal_in);

    IMABC_one_node my_imabc = IMABC_one_node(size_theta,
                           target,
                           target_sd,
                           alpha_start,
                           alpha_goal,
                           N0, Nc, Ngoal, B, LIM1, LIM2, LIM3);
    my_imabc.IMABC_main();
    Rcout << "ESS = " << my_imabc.ESS << "\n";
    print_point_vec(my_imabc.ret);
    write_csv_vpoints(my_imabc.ret, "finalparamsonenode.csv", size_theta, target.size());
}

//' @export
// [[Rcpp::export]]
std::vector<double> check_sim_one_node(NumericVector x){
    //Sim s(80, 2, 5.0, 60, 30, 5, 5);
    Sim_one_node s(x(0), x(1), x(2), x(3), x(4), x(5), x(6));
    return s.read_csv_and_schedule_cancer();
}



