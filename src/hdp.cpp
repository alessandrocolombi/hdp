#include "hdp.h"

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
HdpSampler::HdpSampler( const std::vector<std::vector<double>> &_data,
                        double _priorMean, double _priorA, double _priorB, double _priorLambda,
                        double _a_gamma, double _b_gamma, double _a_alpha, double _b_alpha,
                        double _alpha_init,double _gamma_init,bool _UpdateConc, bool _precompute_Stirling): data(_data), 
                            priorMean(_priorMean), priorA(_priorA), priorB(_priorB), priorLambda(_priorLambda), 
                            a_gamma(_a_gamma), b_gamma(_b_gamma), a_alpha(_a_alpha), b_alpha(_b_alpha),
                            alpha_init(_alpha_init), gamma_init(_gamma_init), UpdateConc(_UpdateConc), precompute_Stirling(_precompute_Stirling)
{
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

    // Compute Stirling numbers
    if(precompute_Stirling){
      Rcpp::Rcout<<"Stirling numbers precalculations ... ";
      lastirling_mat = lastirlings1(numdata); // sub-optimal, it is enough to compute lastirlings1( max_j(n_j) )
      Rcpp::Rcout<<"done!"<<std::endl;  
    }

}


void HdpSampler::init() {

    // Sampling 
    //sample::GSL_RNG engine(seed);  //engine with seed
    sample::rnorm  Norm;           //create normal sampler object
    sample::runif  Unif;           //create uniform sampler object
    sample::sample_index sample_index;  //create categorical sampler object

    //priorMean = 0.0;
    //priorA = 2.0;
    //priorB = 2.0;
    //priorLambda = 0.5;
    alpha = alpha_init;//0.001; // initial value
    gamma = gamma_init;//0.001; // initial value

    int HDP_CLUSTER_RATIO = 5;
    int numClus = numdata / (HDP_CLUSTER_RATIO + numGroups);
    numComponents = numClus;
    TotalNumTables = -1.0; // set in impossible first value, it is never actually used
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
            if( probas.size() == 0 )
                throw std::runtime_error("Error in hdp::init(), probas is empty ");
            if( std::abs( probas.sum()-1 ) > 1e-8 ){
                Rcpp::Rcout<<"probas:"<<std::endl<<probas<<std::endl;
                throw std::runtime_error("Error in hdp::init(), probas do not sum up to 1 ");
            }
            for (int hh=0; hh < probas.size(); hh++) {
                if( probas(hh) < -1e-6 ){
                    Rcpp::Rcout<<"probas("<<hh<<") = "<<probas(hh)<<std::endl;
                    throw std::runtime_error("Error in hdp::init(), probas has negative component ");
                }
            }
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
    //sample::GSL_RNG engine(seed);  //engine with seed
    sample::rnorm  Norm;           //create normal sampler object
    sample::rgamma Gamma;           //create gamma sampler object

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

    // Sampling 
    //sample::GSL_RNG engine(seed);  //engine with seed
    sample::rnorm  Norm;           //create normal sampler object
    sample::sample_index sample_index;  //create categorical sampler object
    sample::rgamma Gamma;           //create gamma sampler object

    double two_pi = 2.0 * M_PI; // 2*\pi

    for (int i = 0; i < numGroups; i++) {
        for (int j=0; j < samplesPerGroup[i]; j++) {
            Eigen::VectorXd logprobas(numComponents + 1);
            int oldAlloc = cluster_allocs[i][j];
            sizes_from_rest[i][oldAlloc] -= 1;
            sizes[oldAlloc] -= 1;
            for (int h=0; h < numComponents; h++) {
                double logproba = std::log(1.0 * sizes_from_rest[i][h] + alpha * betas(h));
                logproba += -0.5*std::log(two_pi*stddevs[h]*stddevs[h]) - 
                             0.5*( (data[i][j] - means[h])*(data[i][j] - means[h]) )/(stddevs[h]*stddevs[h]); //stan::math::normal_lpdf(data[i][j], means[h], stddevs[h]); 
                logprobas[h] = logproba;
            }
            logprobas[numComponents] = marginalLogLikeNormalGamma(data[i][j], priorMean, priorA, priorB, priorLambda); 
            double beta = betas[numComponents];
            logprobas[numComponents] += std::log(alpha * beta);

            Eigen::VectorXd probas = logprobas.array().exp() + 1e-6;
            probas /= probas.sum();

            if( probas.size() == 0 )
                throw std::runtime_error("Error in hdp::sampleAllocations(), probas is empty ");
            if( std::abs( probas.sum()-1 ) > 1e-8 ){
                Rcpp::Rcout<<"probas:"<<std::endl<<probas<<std::endl;
                throw std::runtime_error("Error in hdp::sampleAllocations(), probas do not sum up to 1 ");
            }
            for (int hh=0; hh < probas.size(); hh++) {
                if( probas(hh) < -1e-6 ){
                    Rcpp::Rcout<<"probas("<<hh<<") = "<<probas(hh)<<std::endl;
                    throw std::runtime_error("Error in hdp::sampleAllocations(), probas has negative component ");
                }
            }

            int newAlloc = sample_index(engine, probas);//stan::math::categorical_rng(probas, rng) - 1;
            cluster_allocs[i][j] = newAlloc;
            if (newAlloc == numComponents) {
                std::vector<double> params = normalGammaUpdate(std::vector<double>{data[i][j]}, priorMean, priorA, priorB,priorLambda); 
                double tau = Gamma(engine, params[1], 1.0/params[2]);//stan::math::gamma_rng(params[1], params[2], rng);
                double sigmaNorm = 1.0 / std::sqrt(tau * params[3]);
                double mu = Norm(engine,params[0], sigmaNorm);//stan::math::normal_rng(params[0], sigmaNorm, rng); 
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
    //sample::GSL_RNG engine(seed);  //engine with seed
    sample::sample_index sample_index;  //create categorical sampler object
    sample::rdirichlet<Eigen::VectorXd>  Dirichlet; //create Dirichlet sampler object

    Eigen::VectorXd numTables = Eigen::VectorXd::Ones(numComponents + 1);
    for (int i = 0; i < numGroups; i++) {
        Eigen::VectorXd curr = Eigen::VectorXd::Zero(numComponents + 1);
        for (int h=0; h < numComponents; h++) {
            int numCustomers = sizes_from_rest[i][h];
            Eigen::VectorXd probas = Eigen::VectorXd::Zero(numCustomers);

            // Stirling number computation
            Eigen::VectorXd log_stirling1_vec;
            if(precompute_Stirling){
                Eigen::VectorXd temp = lastirling_mat.row(numCustomers);
                Eigen::VectorXd temp2 = temp.head(numCustomers+1);
                log_stirling1_vec = temp2.tail(numCustomers);
            }
            else
                log_stirling1_vec = lastirling1(numCustomers);

            if(log_stirling1_vec.size() < numCustomers )
                throw std::runtime_error("Error in sampleLatent, size of log_stirling1_vec ");

            for (int m=0; m < numCustomers; m++) {
                if(std::isnan(log_stirling1_vec(m))){
                    Rcpp::Rcout<<"log_stirling1_vec:"<<std::endl<<log_stirling1_vec<<std::endl;
                    throw std::runtime_error("Error in hdp.cpp, get a nan in log_stirling1_vec ");
                }

                /*
                // this is unstable, get inf and 0
                double s = std::exp( log_stirling1_vec(m) ); //double s = 1.0 * stirling_first(numCustomers, m+1); 
                double gammas = std::exp(  std::lgamma(alpha * betas(h)) - std::lgamma(alpha * betas(h) + numCustomers)  );
                probas(m) = s * gammas * std::pow(alpha * betas(h), m+1);
                */
                // it is correct to take the m-th entry because lastirling1 fucntion automatically erase the first value |s(n,0)=0|
                double log_s =  log_stirling1_vec(m); //double s = 1.0 * stirling_first(numCustomers, m+1); 
                double log_gammas = std::lgamma(alpha * betas(h)) - std::lgamma(alpha * betas(h) + numCustomers);
                probas(m) = std::exp( log_s + log_gammas ) * std::pow(alpha * betas(h), m+1);
                if(std::isnan(probas(m))){
                    Rcpp::Rcout<<"log_s = "<<log_s<<std::endl;
                    //Rcpp::Rcout<<"log_stirling1_vec:"<<std::endl<<log_stirling1_vec<<std::endl;
                    Rcpp::Rcout<<"nan m is m = "<<m<<std::endl;
                    Rcpp::Rcout<<"alpha = "<<alpha<<std::endl;
                    Rcpp::Rcout<<"log_gammas = "<<log_gammas<<std::endl;
                    Rcpp::Rcout<<"betas:"<<std::endl<<betas<<std::endl;
                    Rcpp::Rcout<<"numCustomers = "<<numCustomers<<std::endl;
                    throw std::runtime_error("Error in hdp.cpp, get a nan in probas ");
                }
            }
            if (probas.sum() > 0) {
                probas /= probas.sum();

                if( probas.size() == 0 )
                    throw std::runtime_error("Error in hdp::sampleLatent(), probas is empty ");
                if( std::abs( probas.sum()-1 ) > 1e-8 ){
                    Rcpp::Rcout<<"probas:"<<std::endl<<probas<<std::endl;
                    throw std::runtime_error("Error in hdp::sampleLatent(), probas do not sum up to 1 ");
                }
                for (int hh=0; hh < probas.size(); hh++) {
                    if( probas(hh) < -1e-6 ){
                        Rcpp::Rcout<<"probas("<<hh<<") = "<<probas(hh)<<std::endl;
                        throw std::runtime_error("Error in hdp::sampleLatent(), probas has negative component ");
                    }
                }

                curr(h) = sample_index(engine, probas);//stan::math::categorical_rng(probas, rng) - 1;
            } else {
                curr(h) = 0;
            }
        }
        numTables += curr;
    }
    numTables(numComponents) = gamma; // yes, the last Dirichlet parameter is gamma
    betas = Dirichlet(engine, numTables);//stan::math::dirichlet_rng(numTables, rng);
    TotalNumTables = numTables.head(numComponents).sum(); // update number of tables
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

void HdpSampler::save()
{
    out_Allocations.push_back(cluster_allocs);
    out_mu.push_back( Rcpp::NumericVector (means.begin(),means.end()) );  //create temporary vector with current values within push_back call. It is as creating a temporary vector and push it back, but more efficient
    out_sigma.push_back(  Rcpp::NumericVector (stddevs.begin(),stddevs.end())  ); //create temporary vector with current values within push_back call. It is as creating a temporary vector and push it back, but more efficient
    out_K.push_back(numComponents);
    out_alpha.push_back(alpha);
    out_gamma.push_back(gamma);


    // Operations for density estimation -- START
    HDP_Traits::MatRow q_it = HDP_Traits::MatRow::Zero(numGroups, numComponents + 1);
    for (int i = 0; i < numGroups; i++) {
        Eigen::VectorXd log_qs(numComponents + 1);
        for (int h=0; h < numComponents; h++) 
            log_qs[h] = std::log(1.0 * sizes_from_rest[i][h] + alpha * betas(h));  
        
        log_qs[numComponents] = std::log(alpha * betas[numComponents]); // for density estimation
        Eigen::VectorXd qs = log_qs.array().exp() + 1e-6; // for density estimation
        qs /= qs.sum();
        //for (int h=0; h < qs.size(); h++) {
            //if( qs(h) < -1e-6 ){
                //Rcpp::Rcout<<"qs("<<h<<") = "<<qs(h)<<std::endl;
                //throw std::runtime_error("Error, qs has negative component ");
            //}
        //}
        q_it.row(i) = qs;
    }
    out_q.push_back(q_it);
    // Operations for density estimation -- END
}

std::vector<double> HdpSampler::normalGammaUpdate(std::vector<double> data,
                                                  double priorMean, double priorA,
                                                  double priorB, double priorLambda) const
{

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


double HdpSampler::marginalLogLikeNormalGamma(double datum, double mean, double a, double b, double lambda) const
{
  
  std::vector<double> params = normalGammaUpdate(std::vector<double>{datum}, mean, a, b, lambda);

  double out = std::lgamma(params[1]) - std::lgamma(a);
  out += a * std::log(b) - params[1] * std::log(params[2]);
  out += 0.5 * (std::log(lambda) - std::log(params[3]));
  out -= M_PI;
  return out;
}


// Rigon - lastirling1
//
// Vector of lenght n such that the element in k-th position is log|s(n,k)|
// e.g., for n=4 we have |s(n,0)| = 0, |s(n,1)| = 6, |s(n,2)|= 11, |s(n,3)| = 6, |s(n,4)| = 1
// the function returns lastirling1(4) = log(6,11,6,1) = (1.791,2.397,1.791,0)
Eigen::VectorXd HdpSampler::lastirling1(int n) const 
{
  if(n == 0){
    Eigen::VectorXd zero = Eigen::VectorXd::Constant(1,0.0); 
    return (zero);
  }

  Eigen::VectorXd LogSk1 = Eigen::VectorXd::Constant(n+1,0.0); 
  Eigen::VectorXd LogSk = Eigen::VectorXd::Constant(n+1,0.0); 

  LogSk1(1) = 0;
  LogSk1(0) = -std::numeric_limits<double>::infinity();
  LogSk(0)  = -std::numeric_limits<double>::infinity();

  for(int i = 2; i <= n; i++){
    for(int j  = 1; j < i; j++){
      LogSk(j) = LogSk1(j) + std::log(i - 1 + std::exp(LogSk1(j-1) - LogSk1(j)));
    }
    LogSk(i)  = 0;
    LogSk1    = LogSk;
  }
  return( LogSk.tail(n) ); //eliminate the first element, i.e, return the last n elements
}

// Rigon - lastirlings1
//
// Matrix of size (n+1)x(n+1) such that log|s(n,k)| is in position (n+1,k+1) (counting from 1)
// e.g., n = 4; Mat = lastirlings1(n); exp(Mat)[n,] = log(6,11,6,1) = (1.791,2.397,1.791,0) (counting from 0) 
Eigen::MatrixXd HdpSampler::lastirlings1(int n){
  double inf = std::numeric_limits<double>::infinity();

  Eigen::MatrixXd LogS = Eigen::MatrixXd::Constant(n+1,n+1,-inf);
  
  // Fill the starting values
  LogS(0,0) = 0;
  LogS(1,1) = 0;
  
  for(int i = 2; i <= n; i++){
    for(int j = 1; j < i; j++){
      LogS(i,j) = LogS(i-1,j) + std::log(i-1 + std::exp(LogS(i-1,j-1) - LogS(i-1,j))); 
    }
    LogS(i,i)  = 0;
  }

  
  return(LogS);
}

void HdpSampler::updateParams(){

    // Sampling 
    //sample::GSL_RNG engine(seed);  //engine with seed
    sample::rbeta Beta;  //create beta sampler object
    sample::rgamma Gamma; //create gamma sampler object
    sample::runif Unif; //create unif sampler object

    double u{0.0};
    double pi{0.0};

    //1) Concentration parameter gamma 
    double eta = Beta(engine, gamma + 1.0, TotalNumTables);
    pi = (a_gamma + (double)numComponents - 1.0) / (a_gamma + (double)numComponents - 1.0 + TotalNumTables * (b_gamma - std::log(eta)) );
    u = Unif(engine);
    if(u < pi){
        //Rcpp::Rcout<<"-----------------------------------------"<<std::endl;
        //Rcpp::Rcout<<"u < pi; pi = "<<pi<<std::endl;
        //Rcpp::Rcout<<"K = "<<(double)numComponents<<std::endl;
        //Rcpp::Rcout<<"m = "<<TotalNumTables<<std::endl;
        //Rcpp::Rcout<<"shape = "<<a_gamma + (double)numComponents<<std::endl;
        //Rcpp::Rcout<<"rate =  "<< b_gamma - std::log(eta) <<std::endl;
        //Rcpp::Rcout<<"mean = "<< (a_gamma + (double)numComponents) / (b_gamma - std::log(eta)) <<std::endl;
        gamma = Gamma( engine, a_gamma + (double)numComponents,  1.0/(b_gamma - std::log(eta)) );
    }
    else{
        //Rcpp::Rcout<<"-----------------------------------------"<<std::endl;
        //Rcpp::Rcout<<"u > pi; pi = "<<pi<<std::endl;
        //Rcpp::Rcout<<"K = "<<(double)numComponents<<std::endl;
        //Rcpp::Rcout<<"m = "<<TotalNumTables<<std::endl;
        //Rcpp::Rcout<<"shape = "<<a_gamma + (double)numComponents - 1.0<<std::endl;
        //Rcpp::Rcout<<"rate =  "<< b_gamma - std::log(eta) <<std::endl;
        //Rcpp::Rcout<<"mean = "<< (a_gamma + (double)numComponents - 1.0) / (b_gamma - std::log(eta)) <<std::endl;
        gamma = Gamma( engine, a_gamma + (double)numComponents - 1.0,  1.0/(b_gamma - std::log(eta)) );
    }

    //2) Concentration parameter alpha
    double sum_log_w = 0.0;
    unsigned int sum_v = 0;

    for(std::size_t j = 0; j < numGroups; j++){
        sum_log_w += std::log(Beta(engine, alpha + 1.0, (double)samplesPerGroup[j]));
        pi = ( (double)samplesPerGroup[j] )/( (double)samplesPerGroup[j] + alpha);
        u = Unif(engine);
        if(u < sum_v){
            sum_v += 1;
        }
    }
    alpha = Gamma(a_alpha + TotalNumTables - (double)sum_v, 1.0/(b_alpha - sum_log_w));    
}








