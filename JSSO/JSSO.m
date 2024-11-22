function  [Z] =  JSSO(Z1,HSI,MSI,R,p,mu,sd,psf)

[M,N,b] = size(MSI);
X = hyperConvert2D(HSI);
Y = hyperConvert2D(MSI);

[D,S,V] = svds(Z1,p);
C = S*V';
dd1 = zeros(size(D));
eta=1e-6;
mu1=1e-3;
K = 6;
patchsize = 7;
overlap = 4;
bparams.block_sz = [patchsize, patchsize];
bparams.overlap_sz = [overlap overlap];
[nr, nc,~] = size(MSI);
bparams.sz = [nr nc];
sz = [nr nc];


step=patchsize-overlap;
sz1=[1:step:sz(1)- bparams.block_sz(1)+1];
sz1=[sz1 sz(1)- bparams.block_sz(1)+1];
sz2=[1:step:sz(2)- bparams.block_sz(2)+1];
sz2=[sz2 sz(2)- bparams.block_sz(2)+1];
bparams.block_num(1)=length(sz1);
bparams.block_num(2)=length(sz2);
predenoised_blocks = ExtractBlocks1(MSI, bparams);
Y2=Unfold(predenoised_blocks,size(predenoised_blocks),4);
if K==1
    aa = ones(num1*num2,1);
else
  [aa ] = fkmeans(Y2,K);
end
dc1 = zeros(size(C));
C1 = C;  
D1 = D;
for tern= 1:2
        D = D1;
        RDTRD =(R*D)'*(R*D);
	    C = pinv(RDTRD+mu*eye(size(RDTRD,1)))*((R*D)'*Y+mu*(C1-dc1));
        C3 = hyperConvert3D(C1,M,N);
        predenoised_blocks2 = ExtractBlocks1(C3, bparams);
        Z_block2=zeros( bparams.block_sz(1), bparams.block_sz(2),p+(tern-1), bparams.block_num(1)* bparams.block_num(2));
        predenoised_blocks2=permute(predenoised_blocks2,[4 3 1 2]);
        for mn=1:max(aa)
            gg=find(aa==mn);
            XES=predenoised_blocks2(gg,:,:,:);
            [a, b, c, d ]=size(XES);
            XES = reshape(XES,[a b c*d]);   
            V22=Log_prox_tnn( XES, eta/2/mu1 );
            V22=reshape(V22,[a b c d]); 
            Z_block2(:,:,:,gg)=permute(V22,[3 4 2 1]);
        end
        V2= JointBlocks2(Z_block2, bparams);
        V2=V2(1:nr,1:nc,:);
        C1=hyperConvert2D(V2);
        C1 = (C1 + mu*(C+dc1))/(1+mu);
        dc1 = dc1-  (C1-C);   
        C3 = hyperConvert3D(C,M,N);
        CB = imfilter(C3,psf ,'same');
        CBS = CB(1:sd:end, 1:sd:end,:);
        CBS = hyperConvert2D(CBS);
        CBSCBST = CBS*CBS';
        D = (X*CBS'+mu*(D1-dd1))*MoorePenrose(CBSCBST+mu*eye(size(CBSCBST,1)));
        D1 = (D+mu*D1);
        D1 = (D1 + mu*(D+dd1))/(1+mu);
        dd1 = dd1-(D1-D);
        dh  = X-D*CBS;
        [d1,~,~] = svds(dh,1);
        Dn = zeros(size(D,1),size(D,2)+1);
        ddd = zeros(size(Dn));
        ddd(:,1:size(D,2)) = dd1;
        dd1 = ddd; 
        Dn(:,1:size(D,2)) = D;   
        Dn(:,size(D,2)+1) = d1; 
        dz = pinv(R)*Y - D*C;
        [~,k2,l] = svds(dz,1);
        c1 =k2*l';
%         c1 = abs(c1); %% run in matlab2014 need   
        Cn = zeros(size(C,1)+1,size(C,2));
        ddc = zeros(size(Cn));
        ddc(1:size(C,1),:) = dc1;
        dc1 = ddc;      
        Cn(1:size(C,1),:) = C;    
        Cn(size(C,1)+1,:) = c1;  
        C1 = Cn;
        D1 =Dn;
end
Z = hyperConvert3D(D*C,M,N);