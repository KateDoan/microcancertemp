#ifndef SIM_ONE_NODE_H
#define SIM_ONE_NODE_H

// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <math.h>
#include <fstream>
using namespace Rcpp;

class Sim_one_node {
public:
    Sim_one_node(double m1, double mt, double sd_death_other,
                 double m_clinical, double s_clinical,
                 double m_surv_cancer, double s_surv_cancer):
    m1(m1), mt(mt), sd_death_other(sd_death_other),
    m_clinical(m_clinical), s_clinical(s_clinical),
    m_surv_cancer(m_surv_cancer), s_surv_cancer(s_surv_cancer) {};

    void initialize_maps(){
        for(auto year : years){
            map_cancer_incidence[year] = 0;
            map_cancer_death[year] = 0;
            map_other_death[year] = 0;
        }
    }

    void print_map(const std::map<int, int> &m){
        for(auto element : m){
            Rcout << element.first << " " << element.second << "\n";
        }
    }

    void print_num_maps(){
        Rcout << "map_cancer_incidence\n";
        print_map(map_cancer_incidence);
        Rcout << "map_cancer_death\n";
        print_map(map_cancer_death);
        Rcout << "map_other_death\n";
        print_map(map_other_death);
    }

    std::vector<double> concat_map(std::map<int, int> &map_cancer_incidence,
                                   std::map<int, int> &map_cancer_death,
                                   std::map<int, int> &map_other_death){

        std::vector<double> num_vec(3*years.size()-6, 0);
        for(unsigned i=1; i<years.size()-1; i++){
            num_vec[i-1] = (map_cancer_incidence[years[i]] - 26388.0)/30386.0;
        }
        for(unsigned i=1; i<years.size()-1; i++){
            num_vec[i+years.size()-3] = (map_cancer_death[years[i]] - 16563.0)/9486.0;
        }
        for(unsigned i=1; i<years.size()-1; i++){
            num_vec[i+2*years.size()-5] = ((double)map_cancer_death[years[i]]/(double)map_other_death[years[i]] - 0.310616)/0.1078237;
        }
        /*
        for(auto elem : num_vec){
            Rcout << elem << " ";
        }
        Rcout << "\n";
        */

        return num_vec;
    }

    int cal_year(double age_in, int year_in, double age_val){
        return int(year_in + age_val - age_in);
    }

    int find_lower_bound(std::vector<int> vec, int val){
        int res = *(std::upper_bound(vec.begin(), vec.end(), val) - 1);
        return res;
    }

    double rgumbel(double mode, double median){
        double mu = mode;
        double beta = (mode - median)/(std::log( std::log(2.0) ));
        double x = mu - beta * std::log(- std::log(R::runif(0, 1)));
        return x;
    }

    double gen_age_death_rgumbel(double age_low, double mode, double median){
        double age_death = rgumbel(mode, median);
        while(age_low > age_death){
            age_death = rgumbel(mode, median);
            checkUserInterrupt();
        }
        return age_death;
    }

    double gen_age_death_rnorm(double age_low, double mean, double sd){
        double age_death = R::rnorm(mean, sd);
        while(age_low > age_death){
            age_death = R::rnorm(mean, sd);
            checkUserInterrupt();
        }
        return age_death;
    }

    double gen_truncnorm(double age_low, double mean, double sd){
        double age_death;
        int i=0;
        while(true){
            age_death = R::rnorm(mean, sd);
            if(age_death >= age_low) break;
            if(age_low > mean && ++i > 5){
                age_death = age_low + R::runif(0, 5);
                break;
            }
            checkUserInterrupt();
        }
        return age_death;
    }

    void schedule_cancer(double age_in, int year_in){
        double age_clinical, age_death_other, age_death_cancer, age_death;
        int year_clinical = 9999, year_death = 9999;
        std::vector<double> res;

        // parameters
        double mean_death_other = (m1 + (year_in - age_in - 1988) * mt);

        age_death_other = gen_truncnorm(age_in, mean_death_other, sd_death_other);

        int i_cancer = 0;
        while(true){
            age_clinical = gen_truncnorm(0, m_clinical, s_clinical);
            age_death_cancer = gen_truncnorm(age_clinical, age_clinical + m_surv_cancer, s_surv_cancer);
            if(age_death_cancer > age_in){
                break;
            }
            if(++i_cancer > 10){
                age_death_cancer = age_in + R::runif(0, 5);
                break;
            }
            checkUserInterrupt();
        }

        if(age_death_cancer < age_death_other){
            age_death = age_death_cancer;
            year_death = cal_year(age_in, year_in, age_death);
            map_cancer_death[find_lower_bound(years, year_death)] += 1;
        } else {
            age_death = age_death_other;
            year_death = cal_year(age_in, year_in, age_death);
            map_other_death[find_lower_bound(years, year_death)] += 1;
        }

        if(age_clinical < age_death){
            year_clinical = cal_year(age_in, year_in, age_clinical);
            map_cancer_incidence[find_lower_bound(years, year_clinical)] += 1;
        }
        /*
        if(year_death < year_in){
        std::cout << "age_in: " << age_in << " "
                  << "age_clinical: " << age_clinical << " "
                  << "age_death_cancer: " << age_death_cancer << " "
                  << "age_death: " << age_death << "\n"
                  << "year_clinical: " << year_clinical << " "
                  << "year_death: " << year_death << "\n";
        }
        */
    }

    std::vector<double> read_csv_and_schedule_cancer(){
        std::string first_line, line, word;
        std::vector<std::string> row;

        Environment env("package:base");
        Function sysfile = env["system.file"];
        CharacterVector vfname1 = sysfile("extdata", "pop_for_cancer.csv", Named("package", "microcancer"));
        std::string fname1 = as<std::string>(vfname1[0]);
        //Rcout << fname1 << "\n";
        CharacterVector vfname2 = sysfile("extdata", "birth_for_cancer.csv", Named("package", "microcancer"));
        std::string fname2 = as<std::string>(vfname2[0]);
        //Rcout << fname2 << "\n";
        CharacterVector vfname3 = sysfile("extdata", "immigrants_for_cancer.csv", Named("package", "microcancertemp"));
        std::string fname3 = as<std::string>(vfname3[0]);
        //Rcout << fname3 << "\n";

        std::ifstream myfile(fname1);

        if (myfile.is_open())
        {
            getline(myfile, first_line);
            while ( getline (myfile,line) )
            {
                row.clear();
                std::stringstream s(line);
                while (std::getline(s, word, ',')) {
                    row.push_back(word);
                }

                for(unsigned r=0; r<5; r++){
                    for(int i=0; i<std::stoi(row[2])/5; i++){
                        schedule_cancer(std::stoi(row[0]) + r, std::stoi(row[1]));
                    }
                }

                //Rcout << "age: " << row[0] << " year: " << row[1] << " num: " << row[2] << "\n";
            }
            myfile.close();
        }

        else Rcout << "Unable to open file 1";

        myfile.close();

        std::ifstream myfile2(fname2);

        if (myfile2.is_open())
        {
            getline(myfile2, first_line);
            while ( getline (myfile2,line) )
            {
                row.clear();
                std::stringstream s(line);
                while (std::getline(s, word, ',')) {
                    row.push_back(word);
                }

                //for(unsigned r=0; r<5; r++){
                for(int i=0; i<std::stoi(row[2]); i++){
                    schedule_cancer(std::stoi(row[0]), std::stoi(row[1]));
                }
                //}

                //Rcout << "age: " << row[0] << " year: " << row[1] << " num: " << row[2] << "\n";
            }
            myfile2.close();
        }

        else Rcout << "Unable to open file 2";

        myfile2.close();

        std::ifstream myfile3(fname3);

        if (myfile3.is_open())
        {
            getline(myfile3, first_line);
            while ( getline (myfile3,line) )
            {
                row.clear();
                std::stringstream s(line);
                while (std::getline(s, word, ',')) {
                    row.push_back(word);
                }

                //for(int r=0; r<5; r++){
                for(int i=0; i<std::stoi(row[2]); i++){
                    schedule_cancer(std::stoi(row[0]), std::stoi(row[1]));
                }
                //}

                //Rcout << "age: " << row[0] << " year: " << row[1] << " num: " << row[2] << "\n";
            }
            myfile3.close();
        }

        else Rcout << "Unable to open file 3";

        myfile3.close();

        //print_num_maps();

        std::vector<double> numvec = concat_map(map_cancer_incidence, map_cancer_death, map_other_death);
        return numvec;
    }


    private:
        std::map<int, int> map_cancer_incidence;
        std::map<int, int> map_cancer_death;
        std::map<int, int> map_other_death;
        std::vector<int> years = {0, 1988, 1993, 1998, 2003, 2008, 2013};

        double m1, mt, sd_death_other, m_clinical, s_clinical, m_surv_cancer, s_surv_cancer;
};

#endif
