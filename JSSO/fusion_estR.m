function [R] = fusion_estR(HS,MS,sf)
lx = 4;
ly = 4;

[rows2,cols2,bands2] = size(HS);


[nl, nc, bands1] = size(MS);
% Blur operator
middlel = round((nl+1)/2);
middlec = round((nc+1)/2);
% Blur matrix
Bm = zeros(nl, nc);
Bm(middlel-ly:middlel+ly, middlec-lx:middlec+lx) = fspecial('average', [ly*2+1 lx*2+1]);
% Circularly center
Bm = ifftshift(Bm);
% Normalize
Bm = Bm/sum(sum(Bm));
% Fourier transform of the filters
FBm = fft2(Bm);
% Blur MSI
% % Ym = im2mat(MS);
% % Ymb = ConvC(Ym, FBm, nl);
% % % Masked version of Ym (to use when estimating the spatial blur)
% % Ymbim = mat2im(Ymb, nl);
% % LR_MS = downsamp_HS(Ymbim, sf, 1);
% % LR_MS2  = im2mat(LR_MS );

[nlh, nch, L] = size(HS);
% Correspondingly scaled blur's support: [ly*2+1 lx*2+1]
lx = round(lx/sf);
ly = round(ly/sf);
% Blur operator
middlelh=round((nlh+1)/2);
middlech=round((nch+1)/2);
% Blur matrix
Bh=zeros(nlh, nch);
Bh(middlelh-ly:middlelh+ly, middlech-lx:middlech+lx) = fspecial('average', [ly*2+1 lx*2+1]);
% Circularly center
Bh = ifftshift(Bh);
% Normalize
Bh = Bh/sum(sum(Bh));
% Fourier transform of the filters
FBh = fft2(Bh);
% Blur Yh
Yh = im2mat(HS);
HS = ConvC(Yh, FBh, nlh);
% Blur MSI
LR_MS1 = downsamp_HS(MS, sf, 1);
Ym = im2mat(LR_MS1 );
R = zeros(bands1,bands2);
LR_MS2 =  ConvC(Ym , FBh, nlh);
LR_MS = mat2im(LR_MS2, nlh);

% for idx_Lm = 1:bands1
%     no_hs_bands =  length(intersection{idx_Lm});
%     col_aux = zeros(1, no_hs_bands-1);
%     col_aux(1) = 1;
%     row_aux = zeros(1, no_hs_bands);
%     row_aux(1) = 1;
%     row_aux(2) = -1;
%     D = toeplitz(col_aux, row_aux);
%     if (sum(diff(contiguous{idx_Lm})~=1))
%         D(diff(contiguous{idx_Lm})~=1,:) = zeros(sum(diff(contiguous{idx_Lm})~=1), no_hs_bands);
%     end
%     DDt = D'*D;
%     to_inv = HS(intersection{idx_Lm}, :)*HS(intersection{idx_Lm}, :)' + lambda_R*DDt;
%     r = to_inv\(HS(intersection{idx_Lm}, :)*LR_MS2(idx_Lm, :)');
%     R(idx_Lm, intersection{idx_Lm}) = r';
% end
% pre_r=R';

A=HS';

H = A'*A;
options = optimset('Algorithm','interior-point-convex','Display','off','MaxIter',500);

    for k = 1:bands1
        b = double(reshape(LR_MS(:,:,k),[],1));
        f = -1/2*A'*b;
        C = -eye(bands2);
        e = zeros(bands2,1);
%         x0=pre_r(:,k);
        x = quadprog(H,f,C,e,[],[],[],[],[],options);
        R(k,:) = reshape(x,1,[]);

    end 

    Wm=R;
    H2=HS;
    rmse1=0;
    for i=1:100
     Wm=Wm.*(LR_MS2*H2')./(Wm*(H2*H2'));
      aux = sum(sum((LR_MS2 - Wm*H2).^2, 1), 2)/(rows2*cols2);
        rmse2 = sqrt(sum(aux, 3));
        if rmse2-rmse1 < 1e-10
            break;
        end
        rmse1=rmse2;
    end
R=Wm;

    
    

