%-----------------------------------------------------------------------
% Routine for plotting Figure 2(b)
% Author:   Juan Marcos Ramirez Rondon
% Date  :   June, 2017
%
% Ramirez, J. M., & Paredes, J. L. (2014). Robust Sparse Signal Recovery
% Based on Weighted Median Operator. IEEE International Conference on
% Acoustic, Speech, and Signal Processing (ICASSP 2014). pp 1050-1054.
%-----------------------------------------------------------------------
close all;
clear all;

gsim_data = load('nmse_vs_snr_econtaminated.dat');
gpar_data = load('nmse_vs_snr_econtaminated_parameters.dat');

plot(gsim_data(:,1), gsim_data(:,2), 'LineWidth',2);
hold on;
plot(gsim_data(:,1), gsim_data(:,3), 'r','LineWidth',2);
axis('tight'); grid on;
legend('WMR-AR','WMR-HT');
xlabel('SNR [dB]');
ylabel('NMSE [dB]');
title(['Number of realizations = ' num2str(gpar_data(3))]);