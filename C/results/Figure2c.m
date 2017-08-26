%-----------------------------------------------------------------------
% Routine for plotting Figure 2(c)
% Author:   Juan Marcos Ramirez Rondon
% Date  :   June, 2017
%
% Ramirez, J. M., & Paredes, J. L. (2014). Robust Sparse Signal Recovery
% Based on Weighted Median Operator. IEEE International Conference on
% Acoustic, Speech, and Signal Processing (ICASSP 2014). pp 1050-1054.
%-----------------------------------------------------------------------
close all;
clear all;

gsim_data = load('rsnr_vs_snr_implevel.dat');
gpar_data = load('rsnr_vs_snr_implevel_parameters.dat');

plot(gsim_data(:,1)*100, gsim_data(:,2), 'LineWidth',2);
hold on;
plot(gsim_data(:,1)*100, gsim_data(:,3), 'r','LineWidth',2);
axis([0 10 0 25]); grid on;
legend('WMR-AR','WMR-HT');
xlabel('Impulsiveness level [%]');
ylabel('NMSE [dB]');
title(['Number of realizations = ' num2str(gpar_data(3))]);