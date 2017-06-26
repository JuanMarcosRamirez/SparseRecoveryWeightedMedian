/** @file rsnr_vs_snr_implevel.c
 *  @author	Juan Marcos Ramirez Rondon (juanra@ula.ve, juanmarcos26@gmail.com)
 *  @date	June, 2017
 *  @version 	1.00
 *  @brief 	C routine that generates the data for building the Figure 2(c).
 *  @see 	Ramirez, J. M., & Paredes, J. L. (2014). Robust Sparse Signal Recovery Based on Weighted Median Operator. IEEE International Conference on Acoustic, Speech, and Signal Processing (ICASSP 2014). pp 1050-1054. Available on: https://github.com/JuanMarcosRamirez/SparseRecoveryWeightedMedian/blob/master/p1050-ramirez.pdf
 *  @see	Paredes, J. L., & Arce, G. R. (2011). Compressive sensing signal reconstruction by weighted median regression estimates. IEEE Transactions on Signal Processing, 59(6), 2585-2601. Available on: https://www.eecis.udel.edu/~arce/Group/Entries/2012/7/2_Compressive_Spectral_Imaging_files/CompSignReconst.pdf
*/
#include <stdlib.h>
#include <stdio.h>
#include<time.h>
#include<math.h>
#include "../headers/functions.h"


int main()
{
	int num_trials;
	/* Insert number of iterations */
	printf("********************************************************\n");
	printf("Routine that generates data for building the Figure 2(c)\n");
	printf("********************************************************\n");
	printf("Insert the number of realizations per Impulsiveness level: ");
	scanf("%d",&num_trials);

	/* Seed of the pseudo-random number generator
	*/
	srand(time(NULL));
	int i,j,k,l;

	/* Impulsiveness level values
	*/
	int num_points = 11;
	double SNR = 12;
	double imp_level[num_points];
	for(i = 0; i < num_points; i++)
		imp_level[i] = (double)  i / 100.00;

    	/* Sparse signal parameters and projection vector parameters
	*/
    	int N = 1024;
    	int M = 256;
    	int sparsity = 25;

	double x[N];
	double z[M];
	double y[M];
	int *random_index = malloc(N * sizeof(int));
	double energy_x;
	double noise_variance;
	
	/* Dictionary parameters
	*/
	double A[N][M];
	double energyAi;
	
	/* Arrays to save the reconstructed sparse signals
	*/
	double x_hat[N];
	double x_bar[N];

	/* Parameters of the weighted median based algorithms
	*/
	int itmax 	= 100;
	double epsilon 	= 0.01;
	double tol 	= 1e-6;
	double beta 	= 0.95;

	/* Variables for saving the simulation results
	*/
	double mse1, mse2;
	double rsnr_trial_arwmr[num_trials];
	double rsnr_trial_wmr[num_trials];
	double rsnr_snr_arwmr[num_points];
	double rsnr_snr_wmr[num_points];
	double elapsed_time;
	int estimated_time;
	double estimated_sec;
	
	clock_t tic, toc;

	printf("*************************************************************************\n");
	printf("Impulsiveness \tRSNR[dB] \tRSNR[dB] \tElapsed time \tRamaining\n");
	printf("\t\tARWMR\t\tWMR\t\t(sec)\t\ttime \n");
	printf("*************************************************************************\n");
	
	for(k = 0; k < num_points; k++)
	{
		tic = clock();
		for(l = 0; l < num_trials; l++)
		{
		
	        // Building the compressive sensing dictionary
	        // The code exploits the row major property of the C language
			// Therefore, the transpose of the mesurement matrix is built
	        for(i = 0; i < N; i++)
	        {
				energyAi = 0.00;
		        for(j=0; j < M; j++)
		        {
			        A[i][j] = random_normal();
					energyAi += pow(A[i][j],2.0);
		        }

		        for(j=0; j < M; j++)
	        	{
			        A[i][j] = A[i][j] / sqrt(energyAi);
		        }
	        }
		

			// Signal Coefficients and the clean signal energy
	        random_index = random_permutation(N);
	        for(i = 0; i < N; i++)
	        {
		        x[i] = 0.00;
	        }
            energy_x = 0.00;
	        for(i = 0; i < sparsity; i++)
	        {
		        x[random_index[i]] = random_normal();
		        energy_x += (1/(double)N) * pow(x[random_index[i]],2.00);
	        }
		


			// Random Projections
	        for(i = 0; i < M; i++)
	        {
		        z[i] = 0.00;
	        }
	        for(i = 0; i < sparsity; i++)
	        {
		        for(j = 0; j < M; j++)
		        {
		            z[j] += x[random_index[i]] * A[random_index[i]][j];
		        }
	        }
	        
	        // Noisy random projections
	        noise_variance = energy_x / pow(10, SNR/10);
	        for(i = 0; i < M; i++)
	        {
		        y[i] = z[i] + random_econtaminated(imp_level[k], noise_variance, 100);
	        }

		arwmr(N, M, x_hat, y, A, itmax, beta, tol, epsilon);
		wmr(N, M, x_bar, y, A, itmax, beta, tol);

		mse1 = 0.00;
		mse2 = 0.00;
	        for(i = 0; i < N; i++)
	        {
		        mse1 += (1/(double)N) * pow(x[i] - x_hat[i], 2.00);
			mse2 += (1/(double)N) * pow(x[i] - x_bar[i], 2.00);				
	        }
		rsnr_trial_arwmr[l] = - 10 * log10(mse1 / energy_x);
		rsnr_trial_wmr[l] = - 10 * log10(mse2 / energy_x);		
	    }
		toc = clock();
		elapsed_time = (double)(toc - tic) / CLOCKS_PER_SEC;
		estimated_time = (int) (elapsed_time * (num_points - (k + 1))) / 60;
		estimated_sec = (((double)(elapsed_time * (num_points - (k + 1))) / 60) - (double)estimated_time) * 60;
		rsnr_snr_arwmr[k] = 0.00;
		rsnr_snr_wmr[k] = 0.00;
		
		for(i = 0; i < num_trials; i++)
		{
			rsnr_snr_arwmr[k] += (1/(double)num_trials) * rsnr_trial_arwmr[i];
			rsnr_snr_wmr[k] += (1/(double)num_trials) * rsnr_trial_wmr[i];
		}
		printf("%lf\t%lf\t%lf\t%lf\t%d min %d sec\n",imp_level[k], rsnr_snr_arwmr[k], rsnr_snr_wmr[k], elapsed_time, estimated_time, (int)estimated_sec);
	}

	FILE *f = fopen("../results/rsnr_vs_snr_implevel.dat", "w");
	if (f == NULL)
	{
    		printf("Error opening file!\n");
    		exit(1);
	}
	for(i = 0; i < num_points; i++)
		fprintf(f,"%lf\t%lf\t%lf\n",imp_level[i], rsnr_snr_arwmr[i], rsnr_snr_wmr[i]);

	fclose(f);

	f = fopen("../results/rsnr_vs_snr_implevel_parameters.dat", "w");
	if (f == NULL)
	{
    		printf("Error opening file!\n");
    		exit(1);
	}
	fprintf(f,"%d\n",N);
	fprintf(f,"%d\n",M);
	fprintf(f,"%d\n",num_trials);
	fprintf(f,"%d\n",itmax);
	fprintf(f,"%lf\n",beta);
	fprintf(f,"%lf\n",tol);
	fprintf(f,"%lf\n",epsilon);
	fclose(f);

	free(random_index);
	return 0;
}
