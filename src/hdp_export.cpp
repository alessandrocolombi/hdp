// [[Rcpp::depends(RcppEigen)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"
#include "hdp.h"


/*
    // Include Mario
    #include <vector>
    #include <iostream>
    #include <fstream>
    #include <stan/math/prim/mat.hpp>
    
    #include "../hdp.hpp"
    #include "../recordio.hpp"
    #include <random>
    #include "univariate_mixture_state.pb.h"
*/

//' HDP Marginal Sampler
//'
//' @export
// [[Rcpp::export]]
Rcpp::List HDPMarginalSampler(  int Niter, int Nburnin,
                                int d, std::vector<int> n_j, Rcpp::List data_list,
                                double priorMean, double priorA, double priorB, 
                                double priorLambda, double a_gamma, double b_gamma,
                                double a_alpha, double b_alpha,
                                double alpha_init,double gamma_init,bool UpdateConc) {

    int n = std::accumulate(n_j.cbegin(),n_j.cend(),0);

    Rcpp::List& obs = data_list;//Rcpp::as<Rcpp::List>(data_list["observations"]);
    std::vector<std::vector<double>> data;
    data.resize(d);
    // Read data
    for(size_t j = 0; j < d; j++){
        Rcpp::List obs_j = obs[j]; // get list of all individualus at level j
        data[j].reserve(n_j[j]);
        for(size_t i = 0; i < n_j[j]; i++){
            std::vector<double> obs_ji = Rcpp::as<std::vector<double>>(obs_j[i]); // get data for individual i in level j
            data[j].push_back(obs_ji[0]);
        }
    }

    HdpSampler sampler(data,priorMean,priorA,priorB,priorLambda, a_gamma, b_gamma, a_alpha, b_alpha,alpha_init,gamma_init,UpdateConc);
    //Rcpp::Rcout<<"Init start"<<std::endl;
    sampler.init();
    sampler.check();
    //Rcpp::Rcout<<"Init done"<<std::endl;

    Rcpp::Rcout<<"HDP Sampler starts"<<std::endl;
    for (int i=0; i < Nburnin; i++) {
        sampler.sample();
        if( (i+1) % 500 == 0 ) {
            Rcpp::Rcout<<"iter = "<<i<<std::endl;
        }
    }
    
    //Rcpp::Rcout<<"Sampling starts"<<std::endl;
    int thin = 1;
    for (int i=0; i < Niter; i++) {
        sampler.sample();
        if ((i+1) % thin == 0) {
            sampler.save(); //chains.push_back(sampler.getStateAsProto());
        }
        if( (i+1) % 500 == 0 ) {
            Rcpp::Rcout<<"iter = "<<Nburnin + i<<std::endl;
        }
    }

    const std::vector<std::vector< std::vector<int>>>& Alloc = sampler.out_Allocations;
    std::vector<unsigned int> fprowvec; //final partition rowvec
    
    // store total number of data (i.e. cardinality{y_ji})
    Rcpp::NumericMatrix fpmatr(Niter, n );  //final partition matrix

    for (unsigned it = 0; it < Niter; it++){
        fprowvec.clear();
        for(unsigned j=0; j<d; j++){
            fprowvec.insert(fprowvec.end(), Alloc[it][j].cbegin(), Alloc[it][j].cend());
        }
    
        fpmatr(it, Rcpp::_) = Rcpp::NumericMatrix( 1, n, fprowvec.begin() ); //fpmatr[it,:] in R notation
    }
            
    //Rcpp::Rcout << "Done" << std::endl;

    return Rcpp::List::create( Rcpp::Named("Partition") = fpmatr,
                                Rcpp::Named("mu") = sampler.out_mu,
                                Rcpp::Named("sigma") = sampler.out_sigma,
                                Rcpp::Named("alpha") = sampler.out_alpha,
                                Rcpp::Named("gamma") = sampler.out_gamma,
                                Rcpp::Named("K") = sampler.out_K,
                                Rcpp::Named("q_pred") = sampler.out_q
                             ); 
}


//' Testing HDP
//'
//' @export
// [[Rcpp::export]]
Rcpp::List TestHDP() {
    Rcpp::Rcout << "Beginning" << std::endl;
    std::mt19937_64 rng;
    int numGroups = 3;
    int numSamples = 10;//100;

    // Sampling 
    sample::GSL_RNG engine(12345678);  //engine with seed
    sample::rnorm  Norm;           //create normal sampler object
    sample::sample_index sample_index;  //create categorical sampler object
    sample::runif  Unif;           //create uniform sampler object


    Rcpp::Rcout << "numSamples: " << numSamples << std::endl;
    std::vector<std::vector<double>> data(numGroups);
    for (int i=0; i < numGroups; i++) {
        data[i].resize(numSamples);
        for (int j=0; j < numSamples; j++) {
            double u = Unif(engine);//stan::math::uniform_rng(0.0, 1.0, rng);
            if(u < 0.5)
                data[i][j] = Norm(engine,-2.5,1.0);//stan::math::normal_rng(-2.5, 1.0, rng);
            else
                data[i][j] = Norm(engine,2.5,1.0);//stan::math::normal_rng(2.5, 1.0, rng);
        }
    }

    Rcpp::Rcout<<"Stampo data:"<<std::endl;
    for(auto __v : data){
        for(auto __vv : __v)
            Rcpp::Rcout<<__vv<<", ";

        Rcpp::Rcout<<std::endl;
    }
    Rcpp::Rcout<<std::endl;

    /*
        std::ofstream outfile;
        outfile.open("data.csv");
        for (int i = 0; i < numGroups; i++) {
            for (int j=0; j < numSamples; j++) {
                outfile << i << "," << data[i][j] << std::endl;
            }
        }
    */

    //std::deque<HdpState> chains;
    HdpSampler sampler(data,0.0,2.0,2.0,0.5,1.0,1.0,1.0,1.0,1.0,1.0,TRUE);
    Rcpp::Rcout<<"Init start"<<std::endl;
    sampler.init();
    sampler.check();
    Rcpp::Rcout<<"Init done"<<std::endl;
    // spSampler.printDebugString();

    Rcpp::Rcout<<"Burnin starts"<<std::endl;
    int Nburnin = 10;
    for (int i=0; i < Nburnin; i++) {
        sampler.sample();
    }
    // spSampler.printDebugString();
    
    Rcpp::Rcout<<"Sampling starts"<<std::endl;
    int Niter = 10;
    int thin = 1;
    for (int i=0; i < Niter; i++) {
        sampler.sample();
        if ((i+1) % thin == 0) {
            sampler.save(); //chains.push_back(sampler.getStateAsProto());
        }
    }

    const std::vector<std::vector< std::vector<int>>>& Alloc = sampler.out_Allocations;
    std::vector<unsigned int> fprowvec; //final partition rowvec
    
    // store total number of data (i.e. cardinality{y_ji})
    unsigned int n_data = numSamples*numGroups;
    Rcpp::NumericMatrix fpmatr(Niter, n_data );  //final partition matrix

    for (unsigned it = 0; it < Niter; it++){
        fprowvec.clear();
        for(unsigned j=0; j<numGroups; j++){
            fprowvec.insert(fprowvec.end(), Alloc[it][j].cbegin(), Alloc[it][j].cend());
        }
    
        fpmatr(it, Rcpp::_) = Rcpp::NumericMatrix( 1, n_data, fprowvec.begin() ); //fpmatr[it,:] in R notation
    }
            
    // spSampler.printDebugString();
    //writeManyToFile(chains, "chains_hdp.dat");
    Rcpp::Rcout << "Done" << std::endl;

    return Rcpp::List::create( Rcpp::Named("Partition") = fpmatr,
                                Rcpp::Named("mu") = sampler.out_mu,
                                Rcpp::Named("sigma") = sampler.out_sigma,
                                Rcpp::Named("alpha") = sampler.out_alpha,
                                Rcpp::Named("gamma") = sampler.out_gamma,
                                Rcpp::Named("K") = sampler.out_K
                             );
    
    /*
    std::deque<HdpState> restored;
    restored = readManyFromFile<HdpState>("chains_now2.dat");

    Rcpp::Rcout << "******** RESTORED ********" << std::endl;
    */
}