% Per survivor processing
% Reduced state sequence estimation
clear all
close all
clc
data_len = 10^3; % length of the data sequence
num_frames = 10^1; % simulation runs
chan_len = 4; % number of channel taps ( DO NOT CHANGE THIS PARAMETER)
fade_var = 1; % fade variance of the channel
decoding_delay = 10; % decoding delay of the Viterbi algorithm

% SNR parameters
SNR_dB = 30; % SNR per bit (in dB)
SNR = 10^(0.1*SNR_dB);
noise_var = 1*(fade_var*chan_len)/(2*SNR);

C_Ber = 0;% bit errors in each frame initialization
tic()
for i1 = 1: num_frames
% source
a = randi([0 1],1,data_len);
% bpsk mapping
bpsk_seq = 1-2*a;
% impulse response of the ISI channel
fade_chan = normrnd(0,sqrt(fade_var),1,chan_len);
% awgn
noise = normrnd(0,sqrt(noise_var),1,data_len+chan_len-1);
% channel output
chan_op = conv(fade_chan,bpsk_seq)+noise;
% ------------------ RECEIVER----------------------------------------------
% steady state part of the received sequence
steady_state = chan_op(chan_len:data_len);

% branch metrics for the Viterbi algorithm
branch_metric1 = zeros(2^(chan_len-1),data_len-chan_len+1);
branch_metric1(1,:) = (steady_state-(fade_chan(1)+fade_chan(2)+fade_chan(3)+fade_chan(4))).^2;
branch_metric1(2,:) = (steady_state-(fade_chan(1)+fade_chan(2)-fade_chan(3)+fade_chan(4))).^2;
branch_metric1(3,:) = (steady_state-(fade_chan(1)-fade_chan(2)+fade_chan(3)+fade_chan(4))).^2;
branch_metric1(4,:) = (steady_state-(fade_chan(1)-fade_chan(2)-fade_chan(3)+fade_chan(4))).^2;
branch_metric1(5,:) = (steady_state-(-fade_chan(1)+fade_chan(2)+fade_chan(3)+fade_chan(4))).^2;
branch_metric1(6,:) = (steady_state-(-fade_chan(1)+fade_chan(2)-fade_chan(3)+fade_chan(4))).^2;
branch_metric1(7,:) = (steady_state-(-fade_chan(1)-fade_chan(2)+fade_chan(3)+fade_chan(4))).^2;
branch_metric1(8,:) = (steady_state-(-fade_chan(1)-fade_chan(2)-fade_chan(3)+fade_chan(4))).^2;

branch_metric2 = zeros(2^(chan_len-1),data_len-chan_len+1);
branch_metric2(1,:) = (steady_state-(fade_chan(1)+fade_chan(2)+fade_chan(3)-fade_chan(4))).^2;
branch_metric2(2,:) = (steady_state-(fade_chan(1)+fade_chan(2)-fade_chan(3)-fade_chan(4))).^2;
branch_metric2(3,:) = (steady_state-(fade_chan(1)-fade_chan(2)+fade_chan(3)-fade_chan(4))).^2;
branch_metric2(4,:) = (steady_state-(fade_chan(1)-fade_chan(2)-fade_chan(3)-fade_chan(4))).^2;
branch_metric2(5,:) = (steady_state-(-fade_chan(1)+fade_chan(2)+fade_chan(3)-fade_chan(4))).^2;
branch_metric2(6,:) = (steady_state-(-fade_chan(1)+fade_chan(2)-fade_chan(3)-fade_chan(4))).^2;
branch_metric2(7,:) = (steady_state-(-fade_chan(1)-fade_chan(2)+fade_chan(3)-fade_chan(4))).^2;
branch_metric2(8,:) = (steady_state-(-fade_chan(1)-fade_chan(2)-fade_chan(3)-fade_chan(4))).^2;

branch_metric = min(branch_metric1,branch_metric2);

%-----------------------------------------------------------------------------------
% Viterbi algorithm
dec_a=Viterbi_algorithm(data_len-chan_len+1,decoding_delay,branch_metric);
% bit error rate (ignoring transient parts)
C_Ber = C_Ber+nnz(a(4:4+data_len-chan_len-decoding_delay)-dec_a );
end
toc()
% Bit error rate
BER = C_Ber/(num_frames*(data_len-chan_len-decoding_delay+1))