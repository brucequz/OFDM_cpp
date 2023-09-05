
%% OFDM Project
clear;
clc;

figure;
load Pe_bit.mat;
load Pe_symbol.mat;
SNR_db = 0:2:20;
semilogy(SNR_db,Pe_bit(1, :), "-o"); hold on;
semilogy(SNR_db,Pe_bit(2, :), "-x"); hold on;
semilogy(SNR_db,Pe_bit(3, :), "-*"); hold off;
title("BER vs. EbN0 (dB)", 'FontSize', 20);
xlabel('EbN0', 'FontSize',15);
ylabel('Bit Error Rate: P(e)', 'Fontsize', 15);
grid on;
lgd = legend("Bit Error Rate - BPSK", "Bit Error Rate - QPSK", "Bit Error Rate - 16QAM");  
lgd.FontSize = 12;