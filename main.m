clc
clear all
%% load Traffic data
a0=importdata('X14.txt');
for i=1:7
    Z(i,:,:)=a0((i-1)*288+1:i*288,1:144);
end
Z=Z/max(Z(:));
N = size(Z);NN = ndims(Z);
p=0.1;
Omega = find(rand(prod(N),1)<p);O=zeros(size(Z));O(Omega)=Z(Omega);
%% STTF
STTF=STTF_internet(Z,Omega);
time_STTF=toc;
NMAE_STTF = NMAE(Z,STTF,Omega);