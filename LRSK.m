function [result,Z,K,obj,iter]= LRSK(H,y,alpha,beta,gamma)
% initialization
nCluster=length(unique(y));
nn=length(y);
Z=eye(nn);
epsilon = 1e-6;
maxIter = 20;

BTB=H;
obj=[];
for iter = 1:maxIter
    Zold=Z;
    
    %Update Z
    %0.5*tr((-2*Z+Z*Z')*A)+alpha*norm2(Z)^2
    Z= (BTB+2*alpha*eye(nn))\(BTB);
    %Z(find(Z<0))=0;
       
    %Update B
    tmp = H-(eye(nn)-2*Z'+Z*Z')/(2*gamma);
    B = solveB(tmp, gamma/beta);
    BTB=B'*B;
    
    obj(iter)=norm(Z-Zold,'fro');        
    if((iter>5)&(norm(Z-Zold,'fro') < norm(Zold,'fro') * epsilon))
        break
    end
end
K = BTB;
Z = Z - diag(diag(Z));
Z = abs(Z);
Z=(Z+Z')/2;
for ij=1:20
    CKSym = BuildAdjacency(thrC(Z,0.7));
    grps = SpectralClustering(CKSym,nCluster);
    grps = bestMap(y,grps);
    res(ij,:) =  ClusteringMeasure(grps,y);
end
result(1,1)=mean(res(:,1));
result(2,1)=mean(res(:,2));
result(3,1)=mean(res(:,3));
result(1,2)=std(res(:,1));
result(2,2)=std(res(:,2));
result(3,2)=std(res(:,3));
end
