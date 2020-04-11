clear all
close all
clc

addpath('utils')

%Compilation
mex -O CXXFLAGS="\$CXXFLAGS -O3 -Wall" CXXOPTIMFLAGS="-O3" ./SphSPS_mex.cpp -outdir ./

% Image loading
I = double(imread('./data/test_img1.jpg')); 
[h,w,z] = size(I);

% Parameters
K = 1000;  % Superpixel number
m = 0.12;  % Compactness parameter (default 0.12)

% Plug here your contour detector C (see readme.txt file)
C = double(imread('./data/test_img1_contour.png')); 
C = C/max(C(:));

tic;
S = SphSPS_mex(uint8(I), K, m, single(C));
% S = SphSPS(I,K,m);                     %without contour prior
toc;

%% Plotting results
B = spherical_sp_borders(S); %Superpixels borders integrating 360° connectivity
figure,
imagesc(I/255.*repmat(~B,[1 1 3])); 
figure, 
subplot(121); spherical_3D_view(I/255,B,'img')
subplot(122); spherical_3D_view(S,B,'labels')


%% Evaluation of spatial regularity with the Generalized-Global Regularity [1]
GGR = ggr_eval(S)

% For evaluation of more superpixel metrics (ASA, EV, BR, F-measure)
% download the toolbox at:
% https://github.com/rgiraud/sp_toolbox
% and have a look to:
% https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html



