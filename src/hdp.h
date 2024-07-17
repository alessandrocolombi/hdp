#ifndef SRC_HDP_HPP
#define SRC_HDP_HPP

/*
   // Includes from Mario
   //#include <random>
   //#include <vector>
   //#include <stan/math/prim/mat.hpp>

   //#include "mcmc_utils.hpp"
   //#include "stirling_first.hpp"
   //#include "univariate_mixture_state.pb.h"
*/

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"

using namespace HDP_Traits;

/*
   NOTAZIONE:
   P_1,...,P_d | P \iid DP( \alpha, P )
                 P \sim DP( \gamma, P_0 )
               P_0 = Normal( \mu; priorMean, \sigma^2/priorLambda ) x InvGamma( \sigma^2; priorA, priorB )
           \alpha  \sim \operatorname{gamma}(a_\alpha, b_\alpha)
           \gamma  \sim \operatorname{gamma}(a_\gamma, b_\gamma)

*/

class HdpSampler {
 protected:
    int numGroups;
    std::vector<int> samplesPerGroup;
    std::vector<std::vector<double>> data;
    int numdata;

    // mixtures
    int numComponents;
    double TotalNumTables; // this is an integer but in the code it comes from the sum of double. I do not change that
    std::vector<double> means;
    std::vector<double> stddevs;
    std::vector<std::vector<int>> cluster_allocs; // d x n_j vector of vector, same as data
    std::vector<std::vector<int>> sizes_from_rest;
    std::vector<int> sizes;

    // HyperParams for NormalGamma
    double priorMean = 0.0; // default, never used
    double priorA = 2.0; // default, never used
    double priorB = 2.0; // default, never used
    double priorLambda = 0.5; // default, never used

    // HDP
    Eigen::VectorXd betas;
    double alpha = 10;//0.001; // initial value
    double gamma = 1;//0.001; // initial value
    double a_gamma = 2.0; // initial value
    double b_gamma = 2.0; // initial value
    double a_alpha = 2.0; // initial value
    double b_alpha = 2.0; // initial value
    double alpha_init = 0.0;
    double gamma_init = 0.0;
    bool UpdateConc = TRUE;

    // Stirling numbers
    Eigen::MatrixXd lastirling_mat;
    bool precompute_Stirling;

    unsigned long seed = 8012020;
    sample::GSL_RNG engine{seed};
    std::mt19937_64 rng{8012020};


 public:
   // Save
   std::vector<std::vector< std::vector<int>>> out_Allocations;
   std::vector< Rcpp::NumericVector > out_mu;
   std::vector< Rcpp::NumericVector > out_sigma;
   std::vector<int> out_K;
   std::vector<double> out_alpha;
   std::vector<double> out_gamma;
   // WEIGHTS FOR DENSITY ESTIMATION
   // q is a vector of matrices. The external vector has length n_iter. Then, q[it] is a matrix of size d x K+1, where
   // K = K[it] is the number of clusters during iteration it. Within each row, elements are the unnormalized cluster
   // probabilities in log scale for each level j, hence
   // log( q_1(n_11,...,n_1K,U_1) ), ... , log( q_K(n_11,...,n_1K,U_1) ), log( q_{K+1}(n_11,...,n_1K,U_1) )
   // log( q_1(n_j1,...,n_jK,U_j) ), ... , log( q_K(n_j1,...,n_jK,U_j) ), log( q_{K+1}(n_j1,...,n_jK,U_j) )
   // log( q_1(n_d1,...,n_dK,U_d) ), ... , log( q_K(n_d1,...,n_dK,U_d) ), log( q_{K+1}(n_d1,...,n_dK,U_d) )
   std::vector< HDP_Traits::MatRow > out_q;

    ~HdpSampler() = default;

    HdpSampler() {}

    HdpSampler(const std::vector<std::vector<double>> &_data,
               double _priorMean, double _priorA, double _priorB, double _priorLambda,
               double _a_gamma, double _b_gamma,
               double _a_alpha, double _b_alpha,
               double _alpha_init,double _gamma_init,
               bool _UpdateConc, bool _precompute_Stirling);

    void init();

    void sample() {
        sampleAtoms();
        sampleAllocations();
        relabel();
        sampleLatent();
        if(UpdateConc)
         updateParams();
        check();
    }

    void sampleAtoms();

    void sampleAllocations();

    void relabel();

    void sampleLatent();

    void updateParams();

    void check();

    void save();

    std::vector<double> normalGammaUpdate(std::vector<double> data,
                                          double priorMean, double priorA,
                                          double priorB, double priorLambda) const;
    double marginalLogLikeNormalGamma(double datum, double mean, double a, double b, double lambda) const;

    Eigen::VectorXd lastirling1(int n) const;
    Eigen::MatrixXd lastirlings1(int n);
};

#endif  // SRC_HDP_HPP
