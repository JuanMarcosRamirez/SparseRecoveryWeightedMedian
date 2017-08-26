#include <stdio.h>
#include <math.h>
#include "../headers/functions.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Random number generator
// https://stackoverflow.com/questions/7034930/how-to-generate-gaussian-pseudo-random-numbers-in-c-for-a-given-mean-and-varianc
	
double drand()   /* uniform distribution, (0..1] */
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}

double random_normal() 
 /* normal distribution, centered on 0, std dev 1 */
{
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

double random_econtaminated(double level, double variance, double variance_ratio) 
 /* econtaminated distribution distribution, centered on 0, level [0-1], variance */
{
	double u = drand();
	if(u < level)
  		return sqrt(variance_ratio * variance) * sqrt(-2*log(drand())) * cos(2*M_PI*drand());
	else
		return sqrt(variance) * sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

// Random permutation
//https://stackoverflow.com/questions/15961119/how-to-create-a-random-permutation-of-an-array
int* random_permutation(int n) 
{
    int* r = malloc(n * sizeof(int));
    // initial range of numbers
    for(int i=0;i<n;++i){
        r[i]=i+1;
    }
    // shuffle
for (int i = n-1; i >= 0; --i){
    //generate a random number [0, n-1]
    int j = rand() % (i+1);

    //swap the last element with element at random index
    int temp = r[i];
    r[i] = r[j];
    r[j] = temp;
}
  return r;
}

double qselect(double *v, int len, int k)
{
#	define SWAP(a, b) { tmp = v[a]; v[a] = v[b]; v[b] = tmp; }
	int i, st, tmp;
 
	for (st = i = 0; i < len - 1; i++) {
		if (v[i] > v[len-1]) continue;
		SWAP(i, st);
		st++;
	}
 
	SWAP(len-1, st);
 
	return k == st	?v[st]
			:st > k	? qselect(v, st, k)
				: qselect(v + st, len - st, k - st);
}

double absolute_value(double x)
{
	if(x < 0)
		return -x;
	else
		return x;
}


double *array;
static int cmp(const void *a, const void *b){
    int ia = *(int *)a;
    int ib = *(int *)b;
    return array[ia] < array[ib] ? -1 : array[ia] > array[ib];
}

double weighted_median(int M, double samples[], double weights[])
{
	int i,j;
	int indexes[M];

	double threshold;
	double sum_weights = 0.00;
	double sum_temp = 0.00;
	double output;

	for(i = 0; i < M; i++)
	{
		sum_weights += absolute_value(weights[i]);
		indexes[i] = i;
	}
	threshold = 0.50 * sum_weights;	
	
	array = samples;
	qsort(indexes, M, sizeof(*indexes), cmp);	

	j = M-1;	
	while(sum_temp < threshold)
	{
		sum_temp += absolute_value(weights[indexes[j]]);
		j--;
	}	
	output = samples[indexes[j+1]];

	return output;
}



void arwmr(int N, int M, double x[], double y[], double A[][M], int itmax, double beta, double tol, double epsilon)
{
	int i,j;

	double y_aux[M];

	double *absAi = malloc(N * sizeof(double));
	double norm_square_input = 0.00;
	double norm_input;
	double max_absAi;	
	double tau;

	double samples[M + 1];
	double weights[M + 1];
	double x_prev;
	double norm_square_residual;
	double norm_residual;
	
	int iteration = 0;
	double nee = 1.00;

	for(i = 0; i < M; i++)
	{
		y_aux[i] = 0.00;
		norm_square_input += pow(y[i], 2.00);
	}
	norm_input = sqrt(norm_square_input);


	// Computing tau[0]
	for(i = 0; i < N; i++)
	{
		x[i] = 0.00;
		absAi[i] = 0.00;
		for(j = 0; j < M; j++)
		{
			absAi[i] += absolute_value(A[i][j]);
		}
	}
	max_absAi = qselect(absAi, N, N-1);
	tau = 0.50 * epsilon * max_absAi;

	while((iteration < itmax)&&(nee > tol))
	{
		for(i = 0; i < N; i++)
		{
			// Obtaining samples and weights for the i-th coefficient
			for(j = 0; j < M; j++)
			{
				samples[j] = (y[j] - y_aux[j] + x[i] * A[i][j])/A[i][j];
				weights[j] = absolute_value(A[i][j]);
			}
			samples[M] = 0.00;
			weights[M] = tau / (epsilon + absolute_value(x[i]));
			x_prev = x[i];
			x[i] = weighted_median(M+1, samples, weights);
			for(j = 0; j < M; j++)
			{
				y_aux[j] = y_aux[j] + (x[i] - x_prev) * A[i][j];
			}
		}
		
		norm_square_residual = 0.00;
		for(i = 0; i < M; i++)
		{
			norm_square_residual += pow(y[i] - y_aux[i], 2.00); 
		}
		norm_residual = sqrt(norm_square_residual);
		nee = norm_residual / norm_input;
		tau = tau * beta;
		iteration++;
	}
	free(absAi);
}


void wmr(int N, int M, double x[], double y[], double A[][M], int itmax, double beta, double tol)
{
	int i,j;
	

	double y_aux[M];
	double norm_square_input = 0.00;
	double norm_input;
	
	double element_ATy;
	double *abs_ATy = malloc(N * sizeof(double));
	double tau;	

	double residual[M];
	double samples[M];
	double weights[M];
	double x_prev;
	double threshold_compare;
	double norm_square_residual;
	double norm_residual;
	
	int iteration = 0;
	double nee = 1.00;		

	for(i = 0; i < M; i++)
	{
		y_aux[i] = 0.00;
		norm_square_input += pow(y[i], 2.00);
	}
	norm_input = sqrt(norm_square_input);


	// Computing tau[0]
	for(i = 0; i < N; i++)
	{
		x[i] = 0.00;
		element_ATy = 0.00;
		for(j = 0; j < M; j++)
		{
			element_ATy += A[i][j] * y [j];
		}
		abs_ATy[i] = absolute_value(element_ATy);
	}	
	tau = qselect(abs_ATy, N, N-1);


	while((iteration < itmax)&&(nee > tol))
	{
		for(i = 0; i < N; i++)
		{
			// Obtaining samples and weights for the i-th coefficient
			for(j = 0; j < M; j++)
			{
				residual[j] = y[j] - y_aux[j] + x[i] * A[i][j];
				samples[j] = residual[j] / A[i][j];
				weights[j] = absolute_value(A[i][j]);
			}
			x_prev = x[i];
			x[i] = weighted_median(M, samples, weights);	
			
			threshold_compare = 0.00;
			for(j = 0; j < M; j++)
			{
				threshold_compare += absolute_value(residual[j]) - absolute_value(residual[j] - x[i] * A[i][j]) ;
			}
			if(threshold_compare < tau)
				x[i] = 0.00;			


			for(j = 0; j < M; j++)
			{
				y_aux[j] = y_aux[j] + (x[i] - x_prev) * A[i][j];
			}
		}
		
		norm_square_residual = 0.00;
		for(i = 0; i < M; i++)
		{
			norm_square_residual += pow(y[i] - y_aux[i], 2.00); 
		}
		norm_residual = sqrt(norm_square_residual);
		nee = norm_residual / norm_input;
		tau = tau * beta;
		iteration++;
	}
	free(abs_ATy);
}
