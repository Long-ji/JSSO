function  [Z] =  fastfusion(HSI2,MSI2,R,p,sd)
[~,n] = size(MSI2);
[D,S,~] = svds(HSI2,p);
L = size(HSI2,1);
Y = MSI2;
Mz=pinv(R)*Y;
[~,~,V]=svds(Mz,p);
Z = D*(sd*S)*V';
rmse1=0;
for i=1:20
    Z=Z.*(R'*Y)./((R'*R)*Z);
    aux = sum(sum((Y - R*Z).^2, 1), 2)/(n);
    rmse2 = sqrt(sum(aux, 3)/L);
    if rmse2-rmse1 == 0
        break;
    end
    rmse1=rmse2;
end








