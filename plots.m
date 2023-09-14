
%% OFDM Project
clear;
clc;

%%
figure;
load rayleigh_Pebit.mat;
load rayleigh_Pesym.mat;
SNR_db = 0:2:20;
semilogy(SNR_db,Pe_bit(1,:), ":o"); hold on;
semilogy(SNR_db,Pe_bit(2,:), ":x"); hold on;
semilogy(SNR_db,Pe_bit(3,:), ":*"); hold on;

load rayleigh_Pebit_estimate.mat;
load rayleigh_Pesym_estimate.mat;
semilogy(SNR_db,Pe_bit(1,:), "-o"); hold on;
semilogy(SNR_db,Pe_bit(2,:), "-x"); hold on;
semilogy(SNR_db,Pe_bit(3,:), "-*"); hold on;

load rayleigh_Pebit_exact.mat;
load rayleigh_Pesym_exact.mat;
semilogy(SNR_db,Pe_bit(1,:), "--o"); hold on;
semilogy(SNR_db,Pe_bit(2,:), "--x"); hold on;
semilogy(SNR_db,Pe_bit(3,:), "--*"); hold on;

load rayleigh_Pebit_MMSE.mat;
load rayleigh_Pesym_MMSE.mat;
semilogy(SNR_db,Pe_bit(1,:), "-.o"); hold on;
semilogy(SNR_db,Pe_bit(2,:), "-.x"); hold on;
semilogy(SNR_db,Pe_bit(3,:), "-.*"); hold on;

grid on;
title("Rayleigh channel", 'FontSize', 20)
xlabel('EbN0 (dB)', 'FontSize', 15);
ylabel('Bit Error Rate, P(e)', 'FontSize', 15);
lgd = legend("No Equalization - BPSK", "No Equalization - QPSK", "No Equalization - 16QAM", ...
             "ML Estimate - BPSK", "ML Estimate - QPSK", "ML Estimate - 16QAM", ...
             "Exact - BPSK", "Exact - QPSK", "Exact - 16QAM", "MMSE Estimate - BPSK", ...
             "MMSE Estimate - QPSK", "MMSE Estimate - 16QAM");
lgd.FontSize = 12;
%%
figure;
load Pe_bit.mat;
load Pe_symbol.mat;
SNR_db = 0:2:20;
semilogy(SNR_db,Pe_bit(1, :), "-o"); hold on;
semilogy(SNR_db,Pe_bit(2, :), "-x"); hold on;
semilogy(SNR_db,Pe_bit(3, :), "-*"); hold off;
title("Determinstic Channel with full equalization", 'FontSize', 20);
xlabel('EbN0', 'FontSize',15);
ylabel('Bit Error Rate: P(e)', 'Fontsize', 15);
grid on;
lgd = legend("Bit Error Rate - BPSK", "Bit Error Rate - QPSK", "Bit Error Rate - 16QAM");  
lgd.FontSize = 12;