#ifndef __GSLWRAPPERS_H__
#define __GSLWRAPPERS_H__

#include "include_headers.h"
#include "recurrent_traits.h"

//GSL
#include <gsl/gsl_rng.h>     //For random number generators
#include <gsl/gsl_randist.h> //For random variates and probability density functions
#include <gsl/gsl_cdf.h> 	 //For cumulative density functions
#include <gsl/gsl_bspline.h> //For spline operations
#include <gsl/gsl_linalg.h> //For cholesky decomposition
#include <gsl/gsl_sf.h>     //For special functions such as factorials

//Load HDP traits
using namespace HDP_Traits;

/*
	Reference for random number generators:    https://www.gnu.org/software/gsl/doc/html/rng.html
	Reference for random number distributions: https://www.gnu.org/software/gsl/doc/html/randist.html#
*/

namespace sample{ //use the sample:: namespace to avoid clashes with R or other packages

	enum class isChol
	{
		Upper, Lower, False
	};

	/*--------------------------------------------------------------
		Random number generator wrapper
	----------------------------------------------------------------*/

	//This class simply wraps in c++ code the construction and desctruction of a gsl_rng object.
	//I had to remove std::random_device because there is a bug when compiling in window (returs always the same value).
	// reference for bug -> https://en.cppreference.com/w/cpp/numeric/random/random_device, https://sourceforge.net/p/mingw-w64/bugs/338/
	class GSL_RNG{
		public:

			//constructor 1: takes one seed. If seed is 0, generates random seed
			GSL_RNG(unsigned int const & _seed){
				if(_seed == 0){
					seed = static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count());
					std::seed_seq seq = {seed}; //seed provived here has to be random. Than std::seed_seq adds entropy becasuse steady_clock is not sufficientyl widespread
					std::vector<unsigned int> seeds(1);
					seq.generate(seeds.begin(), seeds.end());
					seed = seeds[0];
					//std::cout<<"seed = "<<seed<<std::endl;
				}
				else{
					seed = _seed;
				}
				gsl_rng_env_setup();
				r = gsl_rng_alloc(gsl_rng_default);
				gsl_rng_set(r,seed);
			}

			//constructor 0: default constructor. It is equvalent to the previous one using seed=0.
			GSL_RNG(){
				gsl_rng_env_setup();
				r = gsl_rng_alloc(gsl_rng_default);
				seed = static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count());
				std::seed_seq seq = {seed}; //seed provived here has to be random. Than std::seed_seq adds entropy becasuse steady_clock is not sufficientyl widespread
				std::vector<unsigned int> seeds(1);
				seq.generate(seeds.begin(), seeds.end());
				seed = seeds[0];
				gsl_rng_set(r,seed);
			}

			//descructor, can not be removed because gsl_rng_alloc() is allocating memory
			~GSL_RNG(){
				gsl_rng_free(r);
			}

			//print some information
			void print_info()const{
				printf ("generator type: %s\n", gsl_rng_name(r));
				std::cout<<"seed = "<<seed<<std::endl;
			}

			//call operator, return the engine ready to generate a number
			gsl_rng* operator()()const{
				return r;
			}

			//set seed, not tested much
			inline void set_seed(unsigned int const & s){
				seed = s;
				gsl_rng_set(r,seed);
			}

			//getter
			inline unsigned int get_seed() const{
				return seed;
			}
		private:
			gsl_rng * r; //the random number generator
			unsigned int seed; //the seed
	};

	/*--------------------------------------------------------------
		Random number distribution wrappers
	----------------------------------------------------------------*/

	//Callable object to draw a sample from sampling from Unif([0,1])
	struct runif
	{
		double operator()(GSL_RNG const & engine)const{
			return gsl_rng_uniform(engine()); //gsl_rng_uniform is a function, nothing has to be de-allocated
		}
		double operator()()const{
			return this->operator()(GSL_RNG ());
			//return runif()(GSL_RNG ())

			/*GSL_RNG () creates a GSL_RNG obj calling the default constructor and than calls it call operator.
			In other words, it generates a random number generator, generates a number and destroys the generator.
			It is equivalent to return gsl_rng_uniform( GSL_RNG ()() ).
			In this second case, one () is for the default constructor, the second () is for the call operator.  */
		}
	};

	//Callable object to draw a sample from Unif({0,...,N-1})
	struct runif_int
	{
	  unsigned int operator()(GSL_RNG const & engine, unsigned int const & N)const{
	    return gsl_rng_uniform_int(engine(), N); //gsl_rng_uniform_int is a function, nothing has to be de-allocated.
	  }
	  unsigned int operator()(unsigned int const & N)const{
	    return this->operator()(GSL_RNG (), N);
	    /*GSL_RNG () creates a GSL_RNG obj calling the default constructor and than calls it call operator.
	     In other words, it generates a random number generator, generates a number and destroys the generator*/
	  }
	};
	

	// Generic discrete distribution sampling. Possibly bugged, check it before using it
	struct discrete{

	  	  //Gets the engine
	  //Discrete(p)
	  //OCIO perchè P è un array e non un vector
	    size_t operator()(GSL_RNG const & engine, const double *P) const{
	    size_t K=sizeof(P);
	    //NEEDS PREPROCESSING FOR THE WEIGHTS
	    gsl_ran_discrete_t g=*gsl_ran_discrete_preproc(K, P);
	    return gsl_ran_discrete(engine(),&g);
	  }

	  //Engine defaulted
	  //Discrete (P)
	    size_t operator()(const double *P) const{
	    size_t K=sizeof(P);
	    //NEEDS PREPROCESSING FOR THE WEIGHTS
	    gsl_ran_discrete_t g=*gsl_ran_discrete_preproc(K, P);
	    return gsl_ran_discrete(GSL_RNG ()(),&g);
	  }
	};


	//Callable object to draw a sample from N(mean,sd).
	// --> NB  it takes the standard deviation as input! <--
	struct rnorm
	{
		//Gets the engine
		//N(mean,sd)
		double operator()(GSL_RNG const & engine, double const & mean, double const & sd)const{
			return gsl_ran_gaussian_ziggurat(engine(),sd) + mean;
		}

		//Gets the engine
		//N(0,1)
		double operator()(GSL_RNG const & engine)const{
			return gsl_ran_gaussian_ziggurat(engine(), 1.0);
		}

		//Engine defaulted
		//N(mean,sd)
		double operator()(double const & mean, double const & sd)const{
			return this->operator()(GSL_RNG (), mean, sd);
		}

		//Engine defaulted
		//N(0,1)
		double operator()()const{
			return gsl_ran_gaussian_ziggurat(GSL_RNG ()(),1.0); //the first () is for the constructor, the second il for the call operator
		}
	};

	//Callable object to draw a sample from Gamma(shape,scale).
	// --> NB  Watch out notation! gsl uses the scale parameter, in Europe we are used to the rate parameter (rate = 1/scale) <--
	struct rgamma{

		//Gets the engine
		//Gamma(shape,scale)
		double operator()(GSL_RNG const & engine, double const & shape, double const & scale)const{
			return gsl_ran_gamma(engine(),shape,scale);
		}

		//Engine defaulted
		//Gamma(shape,scale)
		double operator()(double const & shape, double const & scale)const{
			return gsl_ran_gamma(GSL_RNG ()(),shape,scale);
		}
	};

	//Callable object to draw a sample from Chi-squared(k).
	// --> NB  k is the degree of freedom. Chi-squared(k) = Gamma(shape = k/2, scale = 2) <--
	struct rchisq{

		//Gets the engine
		//Chi-squared(k)
		double operator()(GSL_RNG const & engine, double const & k)const{
			return gsl_ran_chisq(engine(),k);
		}

		//Engine defaulted
		//Chi-squared(k)
		double operator()(double const & k)const{
			return gsl_ran_chisq(GSL_RNG ()(), k);
		}
	};

	//Callable object to draw a sample from Poisson(lambda).
	struct rpoisson{

		//Gets the engine
		//Poisson(lambda)
		double operator()(GSL_RNG const & engine, double const & lambda)const{
			return gsl_ran_poisson(engine(),lambda);
		}

		//Engine defaulted
		//Poisson(lambda)
		double operator()(double const & lambda)const{
			return gsl_ran_poisson(GSL_RNG ()(), lambda);
		}
	};

	//Callable object to draw a sample from non-central student-t distribution(n0, mu0, gamma0). The usual t-distribution(n0) can be recovered with mu0=0 and gamma0=1
	struct rnct{

		//Gets the engine
		//non-central student-t distribution(n0, mu0, gamma0)
		double operator()(GSL_RNG const & engine, double const & n0, double const & mu0, const double& gamma0)const{
			if (gamma0 <= 0.0)
				throw std::runtime_error("Error, gamma0 must be strictly positive. ");
			return (gamma0 * gsl_ran_tdist(engine(), n0) + mu0);
		}

		//Engine defaulted
		//non-central student-t distribution(n0, mu0, gamma0)
		double operator()(double const & n0, double const & mu0, const double& gamma0)const{
			return (gamma0 * gsl_ran_tdist(GSL_RNG ()(), n0) + mu0);
		}
	};

	//Callable object to draw a sample from Dirichlet(alpha[0],...,alpha[K-1]).
	//Both input and output Type are template parameters. Return type may be different from Input Type
	// --> NB: Keep record of tested types! It is very challenging to check the type in the code. static_assert does help but undesired behaviours are difficult to be predicted <--
	template<typename RetType = VecCol>
	struct rdirichlet{

		//Default constructor, used to check that Return
		rdirichlet(){
			//Check for Return Type
			static_assert( std::is_same_v<RetType, VecCol> ||
						   std::is_same_v<RetType, VecRow> ||
						   std::is_same_v<RetType, std::vector<double> >  ,
						  "______ ERROR, invalid Return Type requested in rdirichlet. Can handle only VecRow, VecCol and std::vector<double> _____");
		}

		//Gets the engine
		//Dirichlet(alpha[0],...,alpha[K-1])
		template<typename Derived>
		RetType operator()(GSL_RNG const & engine, Eigen::MatrixBase<Derived> const & alpha)const{

			//Check the type of alpha. This makes sure that only eigen vectors are passed or objects that are wrapped into eigen vector by Eigen::Map.
			static_assert( std::is_same_v<Derived, VecCol> ||
						   std::is_same_v<Derived, VecRow> ||
						   std::is_same_v<Derived, Eigen::Map<VecCol>> ||
						   std::is_same_v<Derived, Eigen::Map<const VecCol>> ||
						   std::is_same_v<Derived, Eigen::Map<VecRow>> ||
						   std::is_same_v<Derived, Eigen::Map<const VecRow>> ,
						  "______ ERROR, invalid Input Type requested in rdirichlet. Can handle only VecRow, VecCol  _____");

			unsigned int dim = alpha.size();
			if(dim < 1){
				throw std::runtime_error("length of alpha (concentration parameter) in rdirichlet has to be positive");
			}

			//Declare return object

			RetType return_obj;
			if constexpr(std::is_same_v<RetType, VecCol> || std::is_same_v<RetType, VecRow> ){

				return_obj = RetType::Zero(dim);				//Eigen obj that has to be returned
			}
			else if constexpr(std::is_same_v<RetType, std::vector<double> >){ //static_assert in the constructor avoids any other possibility

				return_obj.resize(dim);
				return_obj = std::vector<double>(dim, 0.0);		//std::vector<double> that has to be returned
			}
			else{
				throw std::runtime_error("______ ERROR, invalid Return Type requested in rdirichlet. Can handle only VecRow, VecCol and std::vector<double> _____");
			}
			gsl_ran_dirichlet(engine(), dim, alpha.derived().data(), return_obj.data());
			//Return
			return(return_obj);
		}

		//Engine defaulted
		//Dirichlet(alpha[0],...,alpha[K-1])
		template<typename Derived>
		RetType operator()(Eigen::MatrixBase<Derived> const & alpha)const{
			return this->operator()(GSL_RNG (), alpha);
		}
	};

	//Callable object to draw a sample from Multinomial(N;K,weights[0],...,weights[K-1]).
	/*  What is a sample from Multinomial(N;K,weights[0],...,weights[K-1]) ?
		It draws a K-dimensional vector (n[0],...,n[K-1]) such that n[0]+...+n[K-1] = N. It is like sampling N times K different colors. Each n[j] counts how may times color j has been selected.
		At each draw, probabilities of the K colors are given by weights.

		Watch out:
		- K is not passed by the user, it is deduced by the length of weights.
		- weights can be normalizied (i.e, they sum up to 1) or not. In the latter case, they are normalized inside gsl_ran_multinomial.
		- it is not checking that weights are non negative and that they are not all 0s.
	*/
	//Both input and output Type are template parameters.
	// --> NB: Keep record of tested types! <--
	template<typename RetType = VecCol>
	struct rmultinomial{

		//Default constructor, used to check that Return
		rmultinomial(){
			//Check for Return Type
			static_assert( std::is_same_v<RetType, VecUnsCol> ||
						   std::is_same_v<RetType, VecUnsRow> ||
						   std::is_same_v<RetType, std::vector<unsigned int> >  ,
						  "______ ERROR, invalid Return Type requested in rmultinomial. Can handle only VecUnsRow, VecUnsCol and std::vector<unsigned int> _____");
		}
		//Gets the engine
		//Multinomial(N;K,weights[0],...,weights[K-1])
		template<typename Derived>
		RetType operator()(GSL_RNG const & engine, unsigned int const & N, Eigen::MatrixBase<Derived> const & weights)const{

			//Check the type of weights. This makes sure that only eigen vectors are passed or objects that are wrapped into eigen vector by Eigen::Map.
			static_assert( std::is_same_v<Derived, VecCol> ||
						   std::is_same_v<Derived, VecRow> ||
						   std::is_same_v<Derived, Eigen::Map<VecCol>> ||
						   std::is_same_v<Derived, Eigen::Map<const VecCol>> ||
						   std::is_same_v<Derived, Eigen::Map<VecRow>> ||
						   std::is_same_v<Derived, Eigen::Map<const VecRow>> ,
						  "______ ERROR, invalid Input Type requested in rmultinomial. Can handle only VecRow, VecCol  _____");

			unsigned int K = weights.size();
			if(K < 1){
				throw std::runtime_error("length of weights in rmultinomial has to be positive");
			}

			//Declare return object

			RetType return_obj;
			if constexpr(std::is_same_v<RetType, VecUnsCol> || std::is_same_v<RetType, VecUnsRow> ){

				return_obj = RetType::Zero(K);				//Eigen obj that has to be returned
			}
			else if constexpr(std::is_same_v<RetType, std::vector<unsigned int> >){ //static_assert in the constructor avoids any other possibility

				return_obj.resize(K);
				return_obj = std::vector<unsigned int>(K, 0.0);		//std::vector<double> that has to be returned
			}
			else{
				throw std::runtime_error("______ ERROR, invalid Return Type requested in rmultinomial. Can handle only VecUnsRow, VecUnsCol and std::vector<unsigned int> _____");
			}
			gsl_ran_multinomial(engine(), K, N, weights.derived().data(), return_obj.data());
			//Return
			return(return_obj);
		}

		//Engine defaulted
		//Multinomial(N;K,weights[0],...,weights[K-1])
		template<typename Derived>
		RetType operator()(unsigned int const & N, Eigen::MatrixBase<Derived> const & weights)const{
			return this->operator()(GSL_RNG (), N,weights);
		}
	};

	// Callable object to draw an index from 0 to M-1 with weights given in weights.
	// Watch out ---> two versions are available.
	// - 1) gets the weights and sets M = weights.size(), i.e the length of the weights decides the number of possible indices
	// - 2) gets the number of indices and sets uniform weights.
	struct sample_index{

		//Gets the engine and the weights
		template<typename Derived>
		unsigned int operator()(GSL_RNG const & engine, Eigen::MatrixBase<Derived> const & weights)const{

			std::vector<unsigned int> n_multinomial_sample = Multinomial(engine, 1, weights); //sample 1 draw form a Multinomial having that weights.

			vector_ui_citerator it( std::find(n_multinomial_sample.cbegin(), n_multinomial_sample.cend(), 1) ); //find the only index that is equal to one.

			if(it == n_multinomial_sample.cend()){ //check
				throw std::runtime_error("Error in sample_index, it is not possible that no level has been selected.");
			}
			else{ //return the position
				return(it - n_multinomial_sample.cbegin());
			}
		}

		//Gets the engine and number of indices
		unsigned int operator()(GSL_RNG const & engine, unsigned int const & M)const{

			if(M <= 0){
				throw std::runtime_error("Error in sample_index, the number of possible indices must be positive");
			}

			VecCol weights = VecCol::Constant(M,1/double(M));
			return this->operator()(engine, weights);
		}

		template<typename Derived>
		unsigned int operator()(Eigen::MatrixBase<Derived> const & weights)const{
			return this->operator()(GSL_RNG (), weights);
		}

		//Gets the engine and number of indices
		unsigned int operator()(unsigned int const & M)const{
			return this->operator()(GSL_RNG (), M);
		}

		private:
			rmultinomial<std::vector<unsigned int>> Multinomial;
	};

	//Multivariate-Normal, Covariance matrix parametrization. Not tested in GDFMM, copied from BGSL
	template<isChol isCholType = isChol::False>
	struct rmvnorm{
		template<typename Derived>
		VecCol operator()(GSL_RNG const & engine, VecCol & mean, Eigen::MatrixBase<Derived>& Cov)
		{
			//Cov has to be symmetric. Not checked for efficiency
			static_assert(isCholType == isChol::False ||
						  isCholType == isChol::Upper ||
						  isCholType == isChol::Lower ,
						  "Error, invalid sample::isChol field inserted. It has to be equal to Upper, Lower or False");
			if(Cov.rows() != Cov.cols())
				throw std::runtime_error("Non squared matrix inserted");
			if(Cov.rows() != mean.size())
				throw std::runtime_error("Matrix and mean do not have compatible size");

			//Declare return object
			/* The simplest way to proceed is to create a gsl_vector for storing the sampled values and then to copy it into an eigen vector. If doing this way the copy cannot be avoided
			   even if Eigen::Map is used. Indeed Eigen map wraps the gsl_vector buffer of data to be an eigen object but when the vector is returned, the gsl_vector has to be freed and the
			   buffer is corrupted, returning the wrong result.
			   What is implemented below is different, first the Eigen object is allocated and then the gsl_vector is defined such that it writes on the same buffer of the Eigen vector. No copies
			   are done before returning.
			   Note that the gsl_vector is not defined by means of gsl_vector_alloc(). Indeed that function allocates some space in memory that is lost when result.data is set to be the same
			   of return_obj, which of course generates a memory leak. This happens because the gsl_vector does not own the buffer and do not deallocate that space by calling gsl_vector_free(). */
			VecCol return_obj(VecCol::Zero(mean.size()));		//Eigen obj that has to be returned
			gsl_vector result;									//gsl_vector where the sampled values will be stored.
			result.size   = mean.size();						//size of the vector
			result.stride = 1;									//how close in memory the elements of result.data are.
			result.owner  = 0;									//result does not own the buffer where data are stored
			result.data   = return_obj.data();					//set data buffer of gsl_vector to the exactly the same of the Eigen vector.
															 	//From now on they share the same buffer, writing on result.data is like writing on return_obj.data() and viceversa

			gsl_matrix *cholMat = gsl_matrix_alloc (Cov.rows(), Cov.rows()); //gsl_ran_multivariate_gaussian requires the Cholesky decompoition and has to be a gsl_matrix
			gsl_vector mu;
			mu.size   = mean.size();
			mu.stride = 1;
			mu.owner  = 0;
			mu.data   = mean.data();
			// Declaration of some Eigen quantities that may be needed. They have to be declared here because they will share their buffer with cholMat
			// and if they go out of scope, the buffer is corrupted.
			MatRow Chol_cov;
			MatRow TCov_row;
			MatCol TCov_col;
			if(Cov.isIdentity()){
				gsl_matrix_set_identity(cholMat);
			}
			else {
				if constexpr( isCholType == isChol::False)
				{
					Chol_cov = Cov.llt().matrixL(); //Use Eigen factorization because Cov is passed by ref and it is not const. If the gsl version is used, Cov is modified!
					cholMat->data = Chol_cov.data();
				}
				else if constexpr(isCholType == isChol::Upper){
					if (Cov.IsRowMajor){
						TCov_row = Cov.transpose(); //Do not use Cov.transposeInPlace() otherwise Cov is modified.
						cholMat->data = TCov_row.derived().data();
					}
					else{
						cholMat->data = Cov.derived().data();
					}
				}
				else if constexpr(isCholType == isChol::Lower){
					if ( !Cov.IsRowMajor ){
						TCov_col = Cov.transpose(); //Do not use Cov.transposeInPlace() otherwise Cov is modified.
						cholMat->data = TCov_col.derived().data();
					}
					else{
						cholMat->data = Cov.derived().data();
					}
				}
			}

			gsl_ran_multivariate_gaussian(engine(), &mu, cholMat, &result);
			//Free and return
			gsl_matrix_free(cholMat);
			return return_obj;
			//Runnig with valgrind_memcheck:
			/*
			==4546== HEAP SUMMARY:
			==4546==     in use at exit: 0 bytes in 0 blocks
			==4546==   total heap usage: 26 allocs, 26 frees, 80,000 bytes allocated
			==4546==
			==4546== All heap blocks were freed -- no leaks are possible
			*/
		}
		template<typename Derived>
		VecCol operator()(VecCol & mean, Eigen::MatrixBase<Derived> & Cov)
		{
			return rmvnorm()(GSL_RNG (), mean, Cov);
		}
	};


	//Follows shape-Scale parametrization. Not tested in GDFMM, copied from BGSL
	template<typename RetType = MatCol, isChol isCholType = isChol::False>
	struct rwish{
		template<typename Derived>
		RetType operator()( GSL_RNG const & engine, double const & b, Eigen::MatrixBase<Derived>& Psi )const
		{
			static_assert(isCholType == isChol::False ||
						  isCholType == isChol::Upper ||
						  isCholType == isChol::Lower ,
						  "Error, invalid sample::isChol field inserted. It has to be equal to utils::isChol::Upper, utils::isChol::Lower or utils::isChol::False");
			if(Psi.rows() != Psi.cols())
				throw std::runtime_error("Non squared matrix inserted");

			RetType return_obj(RetType::Zero(Psi.rows(), Psi.cols()));		//Eigen obj that has to be returned
			gsl_matrix result;												//gsl_matrix where the sampled values will be stored.
			result.size1   = Psi.rows();									//row of the matrix
			result.size2   = Psi.cols();									//cols of the matrix
			result.tda 	  = Psi.rows();										//it is not a submatrix, so this parameter is equal to the number of rows.
			result.owner  = 0;												//result does not own the buffer where data are stored
			result.data   = return_obj.data();								//set data buffer of gsl_matrix to be exactly the same of the Eigen matrix.
															 				//From now on they share the same buffer, writing on result.data is like writing on return_obj.data() and viceversa
			//Declaration of other gsl objects
			gsl_matrix *cholMat = gsl_matrix_calloc(Psi.rows(), Psi.rows());
			gsl_matrix *work    = gsl_matrix_calloc(Psi.rows(), Psi.rows());

			// Declaration of some Eigen quantities that may be needed. They have to be declared here because they will share their buffer with cholMat
			// and if they go out of scope, the buffer is corrupted.
			MatRow Chol_psi;
			MatRow Tpsi_row;
			MatCol Tpsi_col;

			if(Psi.isIdentity()){
				gsl_matrix_set_identity(cholMat);
			}
			else {
				if constexpr( isCholType == isChol::False)
				{
					Chol_psi = Psi.llt().matrixL(); //Use Eigen factorization because Cov is passed by ref and it is not const. If the gsl version is used, Cov is modified!
					cholMat->data = Chol_psi.data();
				}
				else if constexpr(isCholType == isChol::Upper)
				{
					if ( Psi.IsRowMajor ){
						Tpsi_row = Psi.transpose(); //Do not use Psi.transposeInPlace() otherwise Psi is modified.
						cholMat->data = Tpsi_row.derived().data();
					}
					else{
						cholMat->data = Psi.derived().data();
					}
				}
				else if constexpr(isCholType == isChol::Lower)
				{
					if ( !Psi.IsRowMajor ){
						Tpsi_col = Psi.transpose(); //Do not use Cov.transposeInPlace() otherwise Cov is modified.
						cholMat->data = Tpsi_col.derived().data();
					}
					else{
						cholMat->data = Psi.derived().data();
					}
				}
			}

			//Sample with GSL
			gsl_ran_wishart(engine(), b+Psi.rows()-1, cholMat, &result, work);

			//Free and return
			gsl_matrix_free(cholMat);
			gsl_matrix_free(work);
			return return_obj;
			//Running with valgrind
			/*
			==5155== HEAP SUMMARY:
			==5155==     in use at exit: 0 bytes in 0 blocks
			==5155==   total heap usage: 40 allocs, 40 frees, 220,224 bytes allocated
			==5155==
			==5155== All heap blocks were freed -- no leaks are possible
			*/
		}
		template<typename Derived>
		RetType operator()(double const & b, Eigen::MatrixBase<Derived> & Psi)const{
			return rwish<RetType, isCholType>()(GSL_RNG (), b,Psi);
		}
	};
}

#endif
