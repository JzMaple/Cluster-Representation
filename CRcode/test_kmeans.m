clear *; clear -global *; clear all; clc; close all;
S = [1 0.4 0.6 0.2; 0.4 1 0.3 0.7; 0.6 0.3 1 0.4; 0.2 0.7 0.4 1];
preIdx = ones(1,3);
ncluster = 2;
id = kMeans(S,ncluster,preIdx);
disp(id);
