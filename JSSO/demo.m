 
clc;
clear all;


%% load data
addpath(genpath(pwd));
%%%% pavia
datasetName = 'Pavia'; % c 'Pavia' or  'Harvard',

if strcmp(datasetName, 'Pavia')
    load('./data/original_rosis.mat'); 
    S = X(1:256,1:256,:);
elseif strcmp(datasetName, 'Harvard')
    load('./data/imgh7.mat');
    S=ref(1:512,1:512,1:31);
else
    error('Unknown dataset name');
end

S=double(S);
S=S/max(S(:));
HR_HSI = S;
[M, N, B] =size(HR_HSI);
HR_HSI2 = hyperConvert2D(HR_HSI); 

%%   Load spectral response function

if strcmp(datasetName, 'Pavia')
    load('./data/L.mat');
    R=L(:,1:B);
elseif strcmp(datasetName, 'Harvard')
    load('./data/P_N_V2.mat');
    R =P(:,1:B);
else
    error('Unknown dataset name');
end


for i=1:size(R,1)
    sum1=sum(R(i,:));
    for j=1:size(R,2)
        R(i,j)=R(i,j)/sum1;
    end
end
%% simulation HSI and MSI

sd =64; %%%%Spatial downsampling
MSI2 = R*HR_HSI2;
MSI = reshape (MSI2', M,N,size(MSI2,1));
psf =  fspecial('gaussian',7,2);
HR_blur = imfilter(HR_HSI,psf ,'same'); 
HSI = HR_blur(1:sd:end, 1:sd:end,:);
[m, n, B]=size(HSI);
HSI2 = (reshape (HSI,m*n,B))';
MSI2 = hyperConvert2D(MSI);


%%  Real data Hyperion HSI and Sentinel-2A MSI 
% load('HSI.mat')
% load('MSI.mat')

% psf =  fspecial('gaussian',7,2);
% [R] = fusion_estR(HSI,MSI,sd);
% MSI2 = hyperConvert2D(MSI);
% HSI2 =hyperConvert2D(HSI);

%% fusion  
t0=clock;

%%% step 1 Fast Fusion
thr_hsi=fastfusion(HSI2,MSI2,R,1,sd);  %% Scenario I

% thr_hsi= svdsfusion(HSI,MSI,R,psf,sd); %% Scenario II


%%% step 2 Spatial-Spectral Joint Optimization

if strcmp(datasetName, 'Pavia')
    p=7; 
    mu=1e-5;
elseif strcmp(datasetName, 'Harvard')
    p=4; 
    mu=1e-4;
else
    p=7;  
    mu=1e-4;
end


zz = JSSO(thr_hsi,HSI,MSI,R,p,mu,sd,psf);

%%% step 3 Optimization via Error Backpropagation
dHR = imfilter(zz,psf ,'same'); 
dH  = HSI - dHR(1:sd:end, 1:sd:end,:);
method = 'bilinear'; 
H_resized= imresize(dH, sd, method);
H_resized = imresize(H_resized,[M,N]);
out = zz  +  H_resized   ;

t=etime(clock,t0)
[psnr,rmse, ergas, sam, uiqi, ssim, DD, CC] = quality_assessment(double(im2uint8(HR_HSI)), double(im2uint8(out)), 0, 1.0/sd);
fprintf('Quality indices hybrid:\nPSNR=%2.3f REMSE=%2.3f ERGAS=%2.3f SAM=%2.3f UIQI=%2.3f SSIM=%2.3f \n', ...
psnr,rmse,ergas,sam,uiqi,ssim)



