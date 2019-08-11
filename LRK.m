function [result,Z,K,obj,iter]= LRK(H,y,alpha,beta,gamma)
% initialization
nCluster=length(unique(y));
max_mu = 10^3;
nn=length(y);
mu = 0.1;
rho = 20;

Z=eye(nn);
epsilon = 1e-6;
maxIter = 20;
Y=zeros(nn);
K = H;

obj=[];
for iter = 1:maxIter
    Zold=Z;
    %Update Z
    %0.5*tr(K-2*K*Z+Z'*K*Z)+alpha*norm2(Z)^2
    Z= (K+2*alpha*eye(nn))\(K);
    
    %Update W
    G=K-Y/mu;  
    [U, S, V] = svd(G, 'econ');  
    diagS = diag(S);
    svp = length(find(diagS > beta/mu));
    diagS = max(0,diagS - beta/mu);  
    if svp < 0.5
        svp = 1;
    end
    W= U(:,1:svp)*diag(diagS(1:svp))*V(:,1:svp)';   
    %W(find(W<0))=0;    
    
    %Update K
    %0.5*tr(K-2*K*Z+Z'*K*Z)+(gamma/2)*norm2(K-H)^2+(mu/2)*norm2(W-K+Y/mu)^2
    K = mu*W+Y+gamma*H-(0.5*eye(nn)-Z'+0.5*(Z*Z'));
    
    Y = Y+mu*(W-K);
    mu  = min(rho*mu,max_mu);
    mu = 1.1*mu;
    obj(iter)=norm(Z-Zold,'fro');
    if((iter>2)&(norm(Z-Zold,'fro') < norm(Zold,'fro') * epsilon))
        break
    end
end

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
