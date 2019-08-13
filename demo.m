clc;
warning off
ds = {'YALE_165n_1024d_15c_zscore_uni'};

for iData = 1
    dataset = ds{iData}
    data_file = fullfile([dataset, '.mat']);
    kernel_file = fullfile([dataset, '_allkernel.mat']);
    
    load(data_file)
    load(kernel_file)
    
    alpha=1;
    beta=3;
    gamma=10;
    k = 12; %the index of candidate kernel
    H = K(:,:,k);  
    [result,Z,KK1,obj,iter]= LRK(H,y,alpha,beta,gamma);
    %fprintf('#%.3f %.3f %.3f %.3f',alpha,beta,gamma);
    fprintf('# %d> MEA-ACC:%.4f   MEA-NMI:%.4f   MEA-Purity:%.4f\n', rank(KK1), result(1,1), result(2,1),result(3,1));
    fprintf('# %.2f> STD-ACC:%.4f   STD-NMI:%.4f   STD-Purity:%.4f\n', sum(svd(Z)), result(1,2), result(2,2),result(3,2));
    
    [result,Z,KK2,obj,iter]= LRSK(H,y,alpha,beta,gamma);
    fprintf('# %d> MEA-ACC:%.4f   MEA-NMI:%.4f   MEA-Purity:%.4f\n', rank(KK2), result(1,1), result(2,1),result(3,1));
    fprintf('# %.2f> STD-ACC:%.4f   STD-NMI:%.4f   STD-Purity:%.4f\n', sum(svd(KK2)), result(1,2), result(2,2),result(3,2));
    
end