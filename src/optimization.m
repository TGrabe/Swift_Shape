clc;clear;

f=@mainf;
[x,fval,exitflag,output]=fminsearch(f,[-25,50]);