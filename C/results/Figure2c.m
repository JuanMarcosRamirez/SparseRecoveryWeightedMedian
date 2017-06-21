close all;
clear all;

gsim_data = load('rsnr_vs_snr_implevel.dat');
gpar_data = load('rsnr_vs_snr_implevel_parameters.dat');

plot(gsim_data(:,1), gsim_data(:,2), 'LineWidth',2);
hold on;
plot(gsim_data(:,1), gsim_data(:,3), 'r','LineWidth',2);
axis([0 0.2 0 25]); grid on;
legend('ARWM','WMHT');
xlabel('SNR [dB]');
ylabel('NMSE [dB]');
title(['Number of realizations = ' num2str(gpar_data(3))]);