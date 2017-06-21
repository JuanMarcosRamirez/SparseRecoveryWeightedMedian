/** @file nmse_vs_snr_gaussian.c
 *  @author	Juan Marcos Ramirez Rondon (juanra@ula.ve, juanmarcos26@gmail.com)
 *  @date	June, 2017
 *  @version 	1.00
 *  @brief 	C routine that generates the data for building the Figure 2(a).
 *  @see Ramırez, J. M., & Paredes, J. L. (2014). Robust Sparse Signal Recovery Based on Weighted Median Operator. IEEE International Conference on Acoustic, Speech, and Signal Processing (ICASSP 2014). pp 1050-1054. Available on: https://github.com/JuanMarcosRamirez/SparseRecoveryWeightedMedian/blob/master/p1050-ramirez.pdf
 *  @see	Paredes, J. L., & Arce, G. R. (2011). Compressive sensing signal reconstruction by weighted median regression estimates. IEEE Transactions on Signal Processing, 59(6), 2585-2601. Available on: https://www.eecis.udel.edu/~arce/Group/Entries/2012/7/2_Compressive_Spectral_Imaging_files/CompSignReconst.pdf
*/
#include <stdlib.h>
#include <stdio.h>
#include<time.h>
#include<math.h>
#include "../headers/functions.h"


int main()
{
	// Seed of the pseudo-random number generators
	srand(time(NULL));

	int i,j,k,l;

	//Noise parameters
	int num_points = 14;
	double SNR[num_points];
	for(i = 0; i < num_points; i++)
		SNR[i] = 2.00 * i;

	int num_trials = 10;

    	// Signal parameters
    	int N = 512;
    	int M = 256;
    	int sparsity = 25;

	double x[N];
	double x_hat[N];
	double x_bar[N];
	double z[M];
	double y[M];
	int *random_index = malloc(N * sizeof(int));
	double energy_x;
	double noise_variance;
	

	double A[N][M];
	double energyAi;
	
	//Algorithm parameters
	int itmax 	= 100;
	double epsilon 	= 0.01;
	double tol 	= 1e-6;
	double beta 	= 0.95;

	double mse1, mse2;
	double nmse_trial_arwmr[num_trials];
	double nmse_trial_wmr[num_trials];

	double nmse_snr_arwmr[num_points];
	double nmse_snr_wmr[num_points];

	double time_iteration1[num_trials], time_iteration2[num_trials];
	double time_per_iteration1[num_points], time_per_iteration2[num_points];
	
	clock_t tic, toc;
	printf("SNR[dB] \tNMSE[dB] \tmean time \tNMSE[dB] \tmean time\n");
	
	for(k = 0; k < num_points; k++)
	{
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
	        noise_variance = energy_x / pow(10, SNR[k]/10);
	        for(i = 0; i < M; i++)
	        {
		        y[i] = z[i] + sqrt(noise_variance) * random_normal();
	        }

		
		tic = clock();
		arwmr(N, M, x_hat, y, A, itmax, beta, tol, epsilon);
		toc = clock();
		time_iteration1[l] = (double)(toc - tic) / CLOCKS_PER_SEC;

		
		tic = clock();
		wmr(N, M, x_bar, y, A, itmax, beta, tol);
		toc = clock();
		time_iteration2[l] = (double)(toc - tic) / CLOCKS_PER_SEC;

		mse1 = 0.00;
		mse2 = 0.00;
	        for(i = 0; i < N; i++)
	        {
		        mse1 += (1/(double)N) * pow(x[i] - x_hat[i], 2.00);
			mse2 += (1/(double)N) * pow(x[i] - x_bar[i], 2.00);				
	        }
		nmse_trial_arwmr[l] = 10 * log10(mse1 / energy_x);
		nmse_trial_wmr[l] = 10 * log10(mse2 / energy_x);
		
	    }
		
		nmse_snr_arwmr[k] = 0.00;
		time_per_iteration1[k] = 0.00;
		nmse_snr_wmr[k] = 0.00;
		time_per_iteration2[k] = 0.00;
		
		for(i = 0; i < num_trials; i++)
		{
			nmse_snr_arwmr[k] += (1/(double)num_trials) * nmse_trial_arwmr[i];
			time_per_iteration1[k] += (1/(double)num_trials) * time_iteration1[i];
			nmse_snr_wmr[k] += (1/(double)num_trials) * nmse_trial_wmr[i];
			time_per_iteration2[k] += (1/(double)num_trials) * time_iteration2[i];
		}
		printf("%lf\t%lf\t%lf\t%lf\t%lf\n",SNR[k], nmse_snr_arwmr[k], time_per_iteration1[k], nmse_snr_wmr[k], time_per_iteration2[k]);
	}

	FILE *f = fopen("../results/nmse_vs_snr_gaussian.dat", "w");
	if (f == NULL)
	{
    		printf("Error opening file!\n");
    		exit(1);
	}
	for(i = 0; i < num_points; i++)
		fprintf(f,"%lf\t%lf\t%lf\n",SNR[i], nmse_snr_arwmr[i], nmse_snr_wmr[i]);

	fclose(f);

	f = fopen("../results/nmse_vs_snr_gaussian_parameters.dat", "w");
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
