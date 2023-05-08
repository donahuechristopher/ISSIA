%%% Demo for the continuum_removed.m script (see enclosed copyright)
%%% Note: data values should be positive; if negative values are stored in
%%% the input vector, it should be translated to the non-negative orthant

%%% Start

clear all
close all

%%% Generate random data
x = 1:10;
v = randperm(10); 

%%% A way to translate data below
%x = 1:1000;
%v = randn(1,1000);
%v = v+abs(min(v))+eps;

%%% Call continuum_removed function
CR = continuum_removed(v,'Sampling',x,'Plots','yes');

%%% End