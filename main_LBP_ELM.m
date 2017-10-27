%     demo for LBP-ELM classification algorithm
%--------------Brief description-------------------------------------------
%
% 
% This demo implements the LBP-ELM hyperspectral image classification [1]
%
%
% More details in:
%
% [1] W. Li, C. Chen, H. Su, Q. Du. Local Binary Patterns
% and Extreme Learning Machine for Hyperspectral Imagery Classification.
% IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 7, pp.
% 3681-3693, Jul. 2015.
%
% contact: liwei089@ieee.org (Wei Li)

clear all; close all; clc
% 实验1 帕维亚数据集 
load Pavia_University_Data.mat; 
Data0 = z./max(z(:));%像素610*340，103个光谱带
[m n d] = size(Data0);
map = gth;

% band selection 频带选择
X = reshape(Data0, m*n, d);
Psi = PCA_Train(X', 10);% PCA降维
X = X*Psi;
Data = reshape(X, m, n, size(Psi,2));


% LBP feature extraction LBP特征提取
fprintf(' ... ... LBP feature extraction begin ... ...\n');
% LBP的（m，r）参数，可修改（8，2）最优
r = 2;  nr = 8; 
mapping = getmapping(nr,'u2'); 
% LPE频带选择数目-11最优
Feature_P = LBP_feature_global(Data, r, nr, mapping, 11, map);
fprintf(' !!!!!!! LBP feature extraction finished !!!!!!\n');


% save Labelled Data and the labels 保存标签
no_class = max(map(:));
Data = []; Labels = [];
d = size(Feature_P, 3);
Data_tmp = reshape(Feature_P, m*n, d);
for i = 1: no_class
    pos = find(map==i);
    Data = [Data; Data_tmp(pos, :)];
    Labels = [Labels, length(pos)];
end
no_class = length(Labels);%9个标签

%印度数据集16个class
%CTrain = [5 143 83 24 48 73 3 48 2 97 246 59 21 127 39 9];

% 帕维亚数据集9个class
CTrain = [ 66 186 21 31 13 50 13 37 9 ];%426

DataTrain = []; DataTest = [];  CTest = [];
a = 0; 
for i = 1: no_class
    Data_tmp = Data((a+1):(Labels(i)+a), :);
    a = Labels(i) + a;
    rand('seed', 2);
    index_i = randperm(Labels(i));
    DataTrain = [DataTrain; Data_tmp(index_i(1:CTrain(i)), :)];
    DataTest = [DataTest; Data_tmp(index_i(CTrain(i)+1:end), :)];
    CTest =  [CTest length(index_i(CTrain(i)+1:end))];
end
Normalize = max(DataTrain(:));
DataTrain = DataTrain./Normalize;
DataTest = DataTest./Normalize;

% input data 
yapp = []; 
for jj = 1: length(CTrain)
   yapp = [yapp; jj * ones(CTrain(jj),1)];%分类
end
xapp = [yapp DataTrain]; % num_sam x (num_dim + 1) 分类
yapp = []; 
for jj = 1: length(CTest)
   yapp = [yapp; jj * ones(CTest(jj),1)];
end

xtest = [yapp DataTest];%分类

%-----------------------------------------------------
%   Learning Parameters  参数
kerneloption = [5]; %内核选项
c = 1024;
% ELM分类器
[TTrain,TTest,TrainAC,accur_ELM,TY,label] = elm_kernel(xapp,xtest,1,c,'RBF_kernel',kerneloption);
disp(accur_ELM);



