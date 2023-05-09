#include "hdp.hpp"

//HdpSampler::HdpSampler(const std::vector<std::vector<double>> &_data): data(_data) {
    //numGroups = data.size();
    //numdata = 0;
//
    //samplesPerGroup.resize(numGroups);
    //cluster_allocs.resize(numGroups);
    //for (int i=0; i < numGroups; i++) {
        //samplesPerGroup[i] = data[i].size();
        //numdata += samplesPerGroup[i];
        //cluster_allocs[i].resize(samplesPerGroup[i]);
    //}
//}

// Constructor for HdpSampler class that takes in a vector of vectors of doubles as input
HdpSampler::HdpSampler(const std::vector<std::vector<double>> &_data): data(_data) {
    // Set number of groups to size of input data vector
    numGroups = data.size();
    // Initialize number of data points to 0
    numdata = 0;

    // Resize samplesPerGroup and cluster_allocs vectors to the number of groups
    samplesPerGroup.resize(numGroups);
    cluster_allocs.resize(numGroups);

    // For each group in the input data
    for (int i=0; i < numGroups; i++) {
        // Set samplesPerGroup[i] to the size of the i-th group
        samplesPerGroup[i] = data[i].size();
        // Increment numdata by the size of the i-th group
        numdata += samplesPerGroup[i];
        // Resize the cluster_allocs[i] vector to the size of the i-th group
        cluster_allocs[i].resize(samplesPerGroup[i]);
    }
}


void HdpSampler::init() {

    // Sampling 
    sample::GSL_RNG engine(seed);  //engine with seed
    sample::rnorm  Norm;           //create normal sampler object
    sample::runif  Unif;           //create uniform sampler object
    sample::sample_index sample_index;  //create categorical sampler object
    // TODO now we fix this, remember to put priors or pass from
    // constructor
    priorMean = 0.0;
    priorA = 2.0;
    priorB = 2.0;
    priorLambda = 0.5;
    alpha = 1.0;
    gamma = 1.0;

    int HDP_CLUSTER_RATIO = 5;
    int numClus = numdata / (HDP_CLUSTER_RATIO + numGroups);
    numComponents = numClus;
    for (int h=0; h < numClus; h++) {
        means.push_back( Norm(engine,priorMean,0.5) );
        double u{Unif(engine)};
        stddevs.push_back( 0.5 + u * (5.0 - 0.5) );
        //means.push_back(stan::math::normal_rng(priorMean, 5, rng));
        //stddevs.push_back(stan::math::uniform_rng(0.5, 5.0, rng));
    }

    Eigen::VectorXd probas = Eigen::VectorXd::Ones(numClus);
    probas /= probas.sum();

    sizes = std::vector<int>(numClus, 0);
    sizes_from_rest.resize(numGroups);
    for (int i=0; i < numGroups; i++) {
        sizes_from_rest[i] = std::vector<int>(numClus, 0);
        for (int j=0; j < samplesPerGroup[i]; j++) {
            //int clus = stan::math::categorical_rng(probas, rng) - 1;
            int clus = sample_index(engine, probas);
            cluster_allocs[i][j] = clus;
            sizes_from_rest[i][clus] += 1;
            sizes[clus] += 1;
        }
    }

    betas = Eigen::VectorXd::Ones(numClus + 1);
    betas /= betas.sum();
    relabel();
    sampleLatent();
}

void HdpSampler::sampleAtoms() {
    // Sampling 
    sample::GSL_RNG engine(seed);  //engine with seed
    sample::rnorm  Norm;           //create normal sampler object
    sample::rnorm  Gamma;           //create gamma sampler object

    std::vector<std::vector<double>> datavec(numComponents); // one vector per component
    for (int h=0; h < numComponents; h++)
        datavec[h].reserve(numdata);

    //#pragma omp parallel for
    for (int i=0; i < numGroups; i++) {
        for (int j=0; j < samplesPerGroup[i]; j++) {
            int comp = cluster_allocs[i][j];
            datavec[comp].push_back(data[i][j]);
        }
    }

    //#pragma omp parallel for
    for (int h=0; h < numComponents; h++) {
        std::vector<double> params = normalGammaUpdate(datavec[h], priorMean, priorA, priorB, priorLambda);
        double tau = Gamma(engine, params[1], 1.0/params[2] );//stan::math::gamma_rng(params[1], params[2], rng);
        double sigmaNorm = 1.0 / std::sqrt(tau * params[3]);
        double mu = Norm(engine,params[0], sigmaNorm );//stan::math::normal_rng(params[0], sigmaNorm, rng);
        means[h] = mu;
        stddevs[h] = 1.0 / std::sqrt(tau);
    }
    
}

void HdpSampler::sampleAllocations() {
    
    for (int i = 0; i < numGroups; i++) {
        for (int j=0; j < samplesPerGroup[i]; j++) {
            Eigen::VectorXd logprobas(numComponents + 1);
            int oldAlloc = cluster_allocs[i][j];
            sizes_from_rest[i][oldAlloc] -= 1;
            sizes[oldAlloc] -= 1;
            for (int h=0; h < numComponents; h++) {
                double logproba = std::log(1.0 * sizes_from_rest[i][h] + alpha * betas(h));
                logproba += stan::math::normal_lpdf(data[i][j], means[h], stddevs[h]); // <-----------------------------------------------------------------------------
                logprobas[h] = logproba;
            }
            logprobas[numComponents] = utils::marginalLogLikeNormalGamma(data[i][j], priorMean, priorA, priorB, priorLambda); // <-----------------------------------------------------------------------------
            double beta = betas[numComponents];
            logprobas[numComponents] += std::log(alpha * beta);
            Eigen::VectorXd probas = logprobas.array().exp() + 1e-6;
            probas /= probas.sum();
            int newAlloc = stan::math::categorical_rng(probas, rng) - 1; // <-----------------------------------------------------------------------------
            cluster_allocs[i][j] = newAlloc;
            if (newAlloc == numComponents) {
                std::vector<double> params = normalGammaUpdate(std::vector<double>{data[i][j]}, priorMean, priorA, priorB,priorLambda); 
                double tau = stan::math::gamma_rng(params[1], params[2], rng); // <-----------------------------------------------------------------------------
                double sigmaNorm = 1.0 / std::sqrt(tau * params[3]);
                double mu = stan::math::normal_rng(params[0], sigmaNorm, rng); // <-----------------------------------------------------------------------------
                means.push_back(mu);
                stddevs.push_back(1.0 / std::sqrt(tau));
                numComponents += 1;
                sizes.push_back(1);
                for (int k = 0; k < numGroups; k++) {
                    int cnt = (int) k==i;
                    sizes_from_rest[k].push_back(cnt);
                }
                sampleLatent();
            } else  {
                sizes_from_rest[i][newAlloc] += 1;
                sizes[newAlloc] += 1;
            }
        }
    }
    
}

void HdpSampler::sampleLatent() {
    // Sampling 
    sample::GSL_RNG engine(seed);  //engine with seed
    sample::sample_index sample_index;  //create categorical sampler object
    sample::rdirichlet<Eigen::VectorXd>  Dirichlet; //create Dirichlet sampler object

    Eigen::VectorXd numTables = Eigen::VectorXd::Ones(numComponents + 1);
    for (int i = 0; i < numGroups; i++) {
        Eigen::VectorXd curr = Eigen::VectorXd::Zero(numComponents + 1);
        for (int h=0; h < numComponents; h++) {
            int numCustomers = sizes_from_rest[i][h];
            Eigen::VectorXd probas = Eigen::VectorXd::Zero(numCustomers);
            for (int m=0; m < numCustomers; m++) {
                double s = 1.0 * stirling_first(numCustomers, m+1); // <-----------------------------------------------------------------------------
                double gammas = std::exp(
                    std::lgamma(alpha * betas(h) -
                    std::lgamma(alpha * betas(h) + numCustomers)));
                probas(m) = s * gammas * std::pow(alpha * betas(h), m+1);
            }
            if (probas.sum() > 0) {
                probas /= probas.sum();
                curr(h) = sample_index(engine, probas);//stan::math::categorical_rng(probas, rng) - 1;
            } else {
                curr(h) = 0;
            }
        }
        numTables += curr;
    }
    numTables(numComponents) = gamma;
    betas = Dir_veccol(engine, numTables);//stan::math::dirichlet_rng(numTables, rng);
}


void HdpSampler::relabel() {
    std::vector<int> toRemove;
    for (int h=0; h < numComponents; h++) {
        if (sizes[h] == 0)
            toRemove.push_back(h);
    }

    // Remove clusters from the state
    for (auto it = toRemove.rbegin(); it != toRemove.rend(); it++) {
        means.erase(means.begin() + *it);
        stddevs.erase(stddevs.begin() + *it);
        sizes.erase(sizes.begin() + *it);
        for (int i=0; i < numGroups; i++)
            sizes_from_rest[i].erase(sizes_from_rest[i].begin() + *it);
    }

    // adjust allocation labels
    for (int i=0; i < numGroups; i++) {
        for (int j=0; j < samplesPerGroup[i]; j++) {
            int curr = cluster_allocs[i][j];
            cluster_allocs[i][j] -= std::count_if( toRemove.begin(), toRemove.end(),
                                                   [this, &curr](int k) {return k < curr; }
                                                 );
        }
    }
    numComponents = means.size();
    sampleLatent();
}

/*
HdpState HdpSampler::getStateAsProto() {
    HdpState state;
    state.set_num_components(numComponents);
    for (int i=0; i < numGroups; i++) {
        HdpState::GroupParams* p;
        p = state.add_groupparams();
        *p->mutable_cluster_size() = {
            sizes_from_rest[i].begin(), sizes_from_rest[i].end() };
        *p->mutable_cluster_allocs() = {
            cluster_allocs[i].begin(), cluster_allocs[i].end()};
    }
    for (int h=0; h < numComponents; h++) {
        UnivariateMixtureAtom* atom;
        atom = state.add_atoms();
        atom->set_mean(means[h]);
        atom->set_stdev(stddevs[h]);
    }
    *state.mutable_betas() = {betas.data(), betas.data() + betas.size()};
    HdpState::HyperParams hypers;
    hypers.set_mu0(priorMean);
    hypers.set_a(priorA);
    hypers.set_b(priorB);
    hypers.set_lamb(priorLambda);
    state.mutable_hyper_params()->CopyFrom(hypers);
    state.set_alpha(alpha);
    state.set_gamma(gamma);
    return state;
}
*/

void HdpSampler::check() {
    assert (means.size() == stddevs.size());
    assert( means.size() == numComponents);
    std::vector<int> sizes_(numComponents, 0);
    for (int i = 0; i < numGroups; i++) {
        assert(sizes_from_rest[i].size() == numComponents);
        for (int h=0; h < numComponents; h++) {
            sizes_[h] += sizes_from_rest[i][h];
        }
    }
    assert(sizes_ == sizes);
    // std::cout << "numComponents: " << numComponents << std::endl;
    for (int i=0; i < numGroups; i++) {
        std::vector<int> sizes_from_alloc(numComponents, 0);
        // std::cout << "Size: "<< sizes_from_alloc.size() << std::endl;
        for (int j=0; j < samplesPerGroup[i]; j++) {
            // std::cout << "Cluster alloc: " << cluster_allocs[i][j] << std::endl;
            sizes_from_alloc[cluster_allocs[i][j]] += 1;
        }
        assert(sizes_from_alloc == sizes_from_rest[i]);
    }
}


std::vector<double> HdpSampler::normalGammaUpdate(std::vector<double> data,
                                                  double priorMean, double priorA,
                                                  double priorB, double priorLambda) {

  double postMean, postA, postB, postLambda;
  int n = data.size();
  if (n == 0) {
    return std::vector<double>{priorMean, priorA, priorB, priorLambda};
  }

  double sum = std::accumulate(std::begin(data), std::end(data), 0.0);
  double ybar = sum / n;
  postMean = (priorLambda * priorMean + sum) / (priorLambda + n);
  postA = 1.0 * priorA + 1.0 * n / 2;

  double ss = 0.0;
  std::for_each(data.begin(), data.end(),
                [&ss, &ybar](double x) { ss += (x - ybar) * (x - ybar); });

  postB = (priorB + 0.5 * ss +
           0.5 * priorLambda / (n + priorLambda) * n * (ybar - priorMean) *
               (ybar - priorMean));

  postLambda = priorLambda + n;

  return std::vector<double>{postMean, postA, postB, postLambda};
}

