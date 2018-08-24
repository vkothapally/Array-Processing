function [frameindex, nFrames] = getFrameInfo(signal, framesize, framehop)
[nSamples, ~] = size(signal);
nFrames = floor(nSamples/framehop)-1;
frameindex = repmat((1:framesize)',1,nFrames) + (0:nFrames-1)*framehop;



