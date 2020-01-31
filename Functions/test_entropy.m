clc
clear all

%FIND PATH
path0=pwd;

%%ADD REQUIRED SUBFLODERS
path2=[path0 '/shared'];
addpath(path2)
path3=[path0 '/knn'];
addpath(path3)
path4=[path0 '/sqdistance'];
addpath(path4)

%Example 3 (Entropy estimation (base-1: usage))
   Y = rand(5,1000);                         %generate the data of interest (d=5, T=1000)
   mult = 1;                                 %multiplicative constant is important
   co = HShannon_kNN_k_initialization(mult); %initialize the entropy (?H?) estimator
                                              %(?Shannon_kNN_k?), including the value of k
   H = HShannon_kNN_k_estimation(Y,co)      %perform entropy estimation