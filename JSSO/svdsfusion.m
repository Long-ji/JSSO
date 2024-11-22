function  [Z] =  svdsfusion(HSI,MSI,R,psf,sd)
[M, N,~] =size(MSI);
MSI2 = hyperConvert2D(MSI);
HSI2 = hyperConvert2D(HSI);
[~,~,V] = svds(MSI2,1);
[Uh,Sh,~] = svds(HSI2,1);
Z = Uh*sd*Sh*V';
Z3d = hyperConvert3D(Z,M,N);
for i=1:3
    dM = MSI2 - R*Z;
    Z_blur = imfilter(Z3d,psf ,'same'); 
    dH3d = HSI - Z_blur(1:sd:end, 1:sd:end,:);
    dH = hyperConvert2D(dH3d);
    [~,~,V0] = svds(dM,1);
    [Uh,Sh,~] = svds(dH,1);
    Z2 = Uh*sd*Sh*V0';
    Z = Z + Z2;
    Z3d = hyperConvert3D(Z,M,N);
end









