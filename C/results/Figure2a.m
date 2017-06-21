close all;
clear all;

gsim_data = load('nmse_vs_snr_gaussian.dat');
gpar_data = load('nmse_vs_snr_gaussian_parameters.dat');

plot(gsim_data(:,1), gsim_data(:,2), 'LineWidth',2);
hold on;
plot(gsim_data(:,1), gsim_data(:,3), 'r','LineWidth',2);
axis('tight'); grid on;
legend('ARWMH','WMHT');
xlabel('SNR [dB]');
ylabel('NMSE [dB]');
title(['Number of realizations = ' num2str(gpar_data(3))]);