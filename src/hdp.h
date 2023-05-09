#ifndef SRC_HDP_HPP
#define SRC_HDP_HPP

//#include <random>
//#include <vector>
//#include <stan/math/prim/mat.hpp>

//#include "mcmc_utils.hpp"
//#include "stirling_first.hpp"
//#include "univariate_mixture_state.pb.h"

class HdpSampler {
 protected:
    int numGroups;
    std::vector<int> samplesPerGroup;
    std::vector<std::vector<double>> data;
    int numdata;

    // mixtures
    int numComponents;
    std::vector<double> means;
    std::vector<double> stddevs;
    std::vector<std::vector<int>> cluster_allocs;
    std::vector<std::vector<int>> sizes_from_rest;
    std::vector<int> sizes;

    // HyperParams for NormalGamma
    double priorMean, priorA, priorB, priorLambda;

    // HDP
    Eigen::VectorXd betas;
    double alpha;
    double gamma;

    unsigned long seed = 8012020;
    std::mt19937_64 rng{8012020};

 public:
    ~HdpSampler() = default;

    HdpSampler() {}

    HdpSampler(const std::vector<std::vector<double>> &_data);

    void init();

    void sample() {
        sampleAtoms();
        sampleAllocations();
        relabel();
        sampleLatent();
        check();
    }

    void sampleAtoms();

    void sampleAllocations();

    void relabel();

    void sampleLatent();

    void check();

    std::vector<double> normalGammaUpdate(std::vector<double> data,
                                          double priorMean, double priorA,
                                          double priorB, double priorLambda);

    //HdpState getStateAsProto();
};

#endif  // SRC_HDP_HPP