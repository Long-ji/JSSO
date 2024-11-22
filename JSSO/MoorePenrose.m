function  [Z] =  MoorePenrose(G)
[U,S,V]=svd(G);
[m,n]=size(S);
y=min(m,n);
for i=1:y
    if S(i,i)~=0
        S(i,i)=1/S(i,i);
    end
end
S2=zeros(n,m);
for i=1:y
      S2(i,i)=S(i,i);
end
Z=V*S2*U';