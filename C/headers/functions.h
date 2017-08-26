/** @file functions.h
 *  @author	Juan Marcos Ramirez Rondon (juanra@ula.ve, juanmarcos26@gmail.com)
 *  @date	June, 2017
 *  @version 	1.0
 *  @brief Function prototypes for implementing robust sparse signal recovery algorithms based on the weighted median operator.
 *  @see 	Ramirez, J. M., & Paredes, J. L. (2014). Robust Sparse Signal Recovery Based on Weighted Median Operator. IEEE International Conference on Acoustic, Speech, and Signal Processing (ICASSP 2014). pp 1050-1054. Available on: https://github.com/JuanMarcosRamirez/SparseRecoveryWeightedMedian/blob/master/p1050-ramirez.pdf
 *  @see	Paredes, J. L., & Arce, G. R. (2011). Compressive sensing signal reconstruction by weighted median regression estimates. IEEE Transactions on Signal Processing, 59(6), 2585-2601. Available on: https://www.eecis.udel.edu/~arce/Group/Entries/2012/7/2_Compressive_Spectral_Imaging_files/CompSignReconst.pdf
 *
*/
#include<stdlib.h>
#include<stdio.h>

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

/** @return A random number that follows a uniform statistical model between 0 and 1.
*/
double drand();

/** @return A random number obeying to a Gaussian distribution with zero mean and unit variance.
*/
double random_normal(); 

/** Returns a random number following a \f$\varepsilon\f$-contaminated noise model. More precisely, the \f$\varepsilon\f$-contaminated distribution consists in a mixed noise model that merges Gaussian random samples and sparse gross errors. In other words, the \f$\varepsilon\f$-contaminated noise model follows the model given by 
\f[ f_{\eta}(\eta) = \varepsilon \mathcal{N}(0,\sigma^2) + (1-\varepsilon)\mathcal{N}(0,\kappa \sigma^2), \f]
where \f$ \varepsilon\f$ is the impulsiveness level for the mixed noise model;\f$ \sigma^2 \f$ is the variance of the Gaussian noise; and \f$\kappa \f$ is the ratio between the variance of the impulsive noise and the variance of the Gaussian noise.
@return A random number following an e-contaminate noise model.
@param level: Impulsiveness level, \f$\varepsilon\f$.
@param variance: the variance of the Gaussian distribution, \f$\sigma^2\f$.
@param variance_ratio: the ratio between the variance of the impulsive error model and the variance of the Gaussian noise, \f$\kappa\f$. 
*/
double random_econtaminated(double level, double variance, double variance_ratio);

/** @return An array containing a random permutation of integer numbers from 1 to n.
@param n: Number of elements of the output array.
*/
int* random_permutation(int n);

/** Returns the k-order statistic of an array using the quickselect algorithm, where the zero-order statistic is the smallest element of the array and the (N-1)-order statistic is the largest element of the array. 
@return The k-order statistic of the input array
@param v: Input array.
@param len: Length of the input array.
@param k: the required k-order statistic in the range [0,len-1].
*/
double qselect(double *v, int len, int k);

/** @return The absolute value of the input number.
@param x: Input number.
*/
double absolute_value(double x);

static int cmp(const void *a, const void *b);

/** Returns the weighted median of a vector of samples that is weighted by a vector of weights.  
@return The weighted median output.
@param M: Length of the vector of samples (as well as the length of the vector of weights).
@param samples: Vector of samples.
@param weights: Vector of weights.
*/
double weighted_median(int M, double samples[], double weights[]);

/** Returns the sparse coefficient vector given by the adaptive regularizer weighted median regression (arwmr) algorithm.
@return The recovered sparse coefficient vector.
@param N: Coefficient vector length.
@param M: Projection vector length.
@param x[]: output coefficient vector.
@param y[]: Projection vector. 
@param A[][M]: Dicctionary.
@param itmax: Maximum number of iterations.
@param beta: Continuation parameter.
@param tol: Tolerance
@param epsilon: epsilon.
*/
void arwmr(int N, int M, double x[], double y[], double A[][M], int itmax, double beta, double tol, double epsilon);

/** Returns the sparse coefficient vector given by the weighted median regression (wmr) algorithm.
@return The recovered sparse coefficient vector.
@param N: Coefficient vector length.
@param M: Projection vector length.
@param x[]: output coefficient vector.
@param y[]: Projection vector. 
@param A[][M]: Dicctionary.
@param itmax: Maximum number of iterations.
@param beta: Continuation parameter.
@param tol: Tolerance
*/
void wmr(int N, int M, double x[], double y[], double A[][M], int itmax, double beta, double tol);

#endif
