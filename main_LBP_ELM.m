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
% ʵ��1 ��ά�����ݼ� 
load Pavia_University_Data.mat; 
Data0 = z./max(z(:));%����610*340��103�����״�
[m n d] = size(Data0);
map = gth;

% band selection Ƶ��ѡ��
X = reshape(Data0, m*n, d);
Psi = PCA_Train(X', 10);% PCA��ά
X = X*Psi;
Data = reshape(X, m, n, size(Psi,2));


% LBP feature extraction LBP������ȡ
fprintf(' ... ... LBP feature extraction begin ... ...\n');
% LBP�ģ�m��r�����������޸ģ�8��2������
r = 2;  nr = 8; 
mapping = getmapping(nr,'u2'); 
% LPEƵ��ѡ����Ŀ-11����
Feature_P = LBP_feature_global(Data, r, nr, mapping, 11, map);
fprintf(' !!!!!!! LBP feature extraction finished !!!!!!\n');


% save Labelled Data and the labels �����ǩ
no_class = max(map(:));
Data = []; Labels = [];
d = size(Feature_P, 3);
Data_tmp = reshape(Feature_P, m*n, d);
for i = 1: no_class
    pos = find(map==i);
    Data = [Data; Data_tmp(pos, :)];
    Labels = [Labels, length(pos)];
end
no_class = length(Labels);%9����ǩ

%ӡ�����ݼ�16��class
%CTrain = [5 143 83 24 48 73 3 48 2 97 246 59 21 127 39 9];

% ��ά�����ݼ�9��class
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
   yapp = [yapp; jj * ones(CTrain(jj),1)];%����
end
xapp = [yapp DataTrain]; % num_sam x (num_dim + 1) ����
yapp = []; 
for jj = 1: length(CTest)
   yapp = [yapp; jj * ones(CTest(jj),1)];
end

xtest = [yapp DataTest];%����

%-----------------------------------------------------
%   Learning Parameters  ����
kerneloption = [5]; %�ں�ѡ��
c = 1024;
% ELM������
[TTrain,TTest,TrainAC,accur_ELM,TY,label] = elm_kernel(xapp,xtest,1,c,'RBF_kernel',kerneloption);
disp(accur_ELM);



