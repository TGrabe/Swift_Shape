function [n] = n_pmma(lambdanm)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

lambda=lambdanm/1000;
a0=2.1778;
a1=6.1209*10^-3;
a2=-1.5004*10^-3;
a3=2.3678*10^-2;
a4=-4.2137*10^-3;
a5=7.3417*10^-4;
a6=-4.5042*10^-5;
n=sqrt(a0+a1*lambda^2+a2*lambda^4+a3*lambda^-2+a4*lambda^-4+a5*lambda^-6+a6*lambda^-8);
end

