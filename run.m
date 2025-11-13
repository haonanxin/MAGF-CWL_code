clear;clc;close all;
addpath('funs');
addpath('data');

dataset_name='BBC';
load([dataset_name,'.mat'])
v = length(X);
c = length(unique(Y));
n = length(Y);

%% Normalization
for i = 1 :length(X)
    X{i} = full((X{i} - mean(X{i}, 2)) ./ repmat(std(X{i}, [], 2), 1, size(X{i}, 2)));
end

%% Parameter Setting of BBC with normalization               ACC = 0.85
L=3;
alpha=0.1;

%% Parameter Setting of Digit4k with normalization        ACC = 0.91
% L=3;
% alpha=0.3;

%% Parameter Setting of mnist4 with normalization     ACC = 0.90
% L=5;
% alpha=0.8;

%% Parameter Setting of MNIST_mv without normalization       ACC = 0.90
% L=3;
% alpha=0.5;


%% Optimization of MAGF-CWL
[F,num_view] = Initialize_F_FINCH(X,c);
[B] = AGC_FINCH(X,L,alpha);
[learned_F,loss] = main_MAGF_CWL(X,L-1,B,c,F);

[~, ind] = max(learned_F, [], 2);
out = ClusteringMeasure_new(Y, ind);


disp(['********************************************']);
disp(['Running MAGF-CWL on ',dataset_name,' to obtain ACC: ', num2str(out.ACC)]);
disp(['********************************************']);