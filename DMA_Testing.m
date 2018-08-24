% PURPOSE: Understanding Differential Microphonne Array
% Author:  Vinay Kothapally
% Version: 1.0       Date:    08/23/2018
% Revised: 08/24/2018  (Signal recorded and then moved)

%=====================%
clc; close all; clear;
%=====================%


addpath('/home/vkk160330/Documents/Matlab Toolboxes/FastISM');

%=============== Room Simulations ===============%
Room.Fs = 16000;               
Room.room = [10  10  10];      
Room.room_center = Room.room/2;

d = 10e-3; c = 343; M = 4;
x = (-1.5:1.5)'*d; resolution = 5; 
Room.mic_pos = [x 0*x 0*x] + Room.room_center;
Room.T60 = 0;                 
Room.abs_weights = ones(1,6); 

signal = 0.01*randn(1*Room.Fs,1); 
theta = (0:resolution:360-resolution)';
AuData  = zeros(length(signal), M, length(theta));
for k = 1:length(theta)
    Room.src_traj = [sind(theta(k)) cosd(theta(k))  0] + Room.room_center;
    fast_ISM_RIR_bank(Room,'fastISM_RIRs.mat');
    AuData(:,:,k) = ISM_AudioData('fastISM_RIRs.mat',signal);
end

%=============== Frame Processing ===============%
framesize = 256; framehop = 128; nFFT = 256;
[frameindex, nFrames] = getFrameInfo(AuData, 256, 128);
win = squeeze(hamming(framesize)); nChannels = size(AuData,2);


%=============== DMA Filter Design ==============%
f = (0:nFFT/2)'*Room.Fs/nFFT; t0 = d/c; alpha = -1;
f(f<80) = 80;
for k = 1:length(f)
    d = [exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(90)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(0)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(30)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(60)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(120)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(150)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(180)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(220)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(250)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(280)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(310)); ...
         exp(-1i*2*pi*f(k)*(0:M-1)*t0*cosd(340))];          
    H(k,:) = pinv(d)*[1;zeros(length(d)-1,1)];
end

%============ DMA Processing (Freq) ==============%
beam = zeros(frameindex(end,end), length(theta));
for kk = 1:length(theta)
for k = 1:nFrames
    frameData = AuData(frameindex(:,k), :).*(win*ones(1,nChannels));
    frameFFT  = fft(frameData, nFFT); frameFFT = frameFFT(1:nFFT/2+1,:);
    OutFFT    = sum(frameFFT.*H,2); OutFFT = [OutFFT; conj(flipud(OutFFT(2:end-1)))];
    beam(frameindex(:,k)) = beam(frameindex(:,k)) + real(ifft(OutFFT));
    
end
end


%============ Beam Pattern (Polar) ==============%
theta = linspace(0, 360, length(beam))';
[T,R] = cart2pol(beam.*cosd(theta), beam.*sind(theta));
polar(T,R);




rmpath('/home/vkk160330/Documents/Matlab Toolboxes/FastISM');
    
    
    



