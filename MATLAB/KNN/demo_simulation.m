%% KNN Fused algorithm using simulation


% Samples and indices
n= 5000;
T = rand(1,n);  u =T;

% randomly choosing Q (quantile)  
Q = rand(1,n);
Q = 0.9*Q + 0.05; 

% generating data X and Y
X1 = repmat(1,1,n);
X2 = normrnd(0,1,[1,n]);
X3 = rand(1,n);
X4 = rand(1,n); X4 = X4 > 0.5;
X5 = rand(1,n); X5 = X5 > 0.5;
X6 = rand(1,n);  
X7 = rand(1,n);  
X8 = normrnd(0,1,[1,n]); 
X9 = normrnd(0,1,[1,n]);
X = [X1;X2;X3;X4;X5;X6;X7;X8;X9]; X=X';   % n by p
errs = normrnd(0,0.01,[1,n]);
Y = repmat(1,1,n)-7*X(:,2)'+(10+2*sin(2*pi*T)).*X(:,3)'+(3+2*Q).*X(:,4)'+5*(Q>0.5).*X(:,5)'+5*(T>0.5).*X(:,6)'+3*(Q+T).*X(:,7)'+errs;


% underlying varying coefficient function
beta=zeros(n,9);
beta(:,1)= 1+ norminv(Q,0,0.01); 
beta(:,2)= repmat(-7,1,n);
beta(:,3)= (10+2*sin(2*pi*T));
beta(:,4)= 3+2*Q;
beta(:,5)= 5*(Q>0.5);
beta(:,6)= 5*(T>0.5);
beta(:,7)= repmat(0,1,n); 
beta(:,7)= 3*Q+3*T;
beta(:,8)= repmat(0,1,n); beta(:,9)=repmat(0,1,n); 


%% Run the knn admm algorithm with fixed lambda 
% K : nearest neighbor number
% lambda: penalty parameter
% eta: step-size parameter of the ADMM algorithm
% tolerance: tolerance parameter determing stopping criterion of the ADMM 
% Here we set K=5, lambda=1, tolerance =0.1, eta (step-size parameter) = 0.5, 

lambda = 1;
[theta2, theta_set2] = VC_qt_knn_admm(X,Q,T,Y,5,lambda,500, 0.1, 0.5);


%% Underlying and estimated coefficient plots
color_limits=[];
for ii = 1:9
subplot(3,3,ii)
scatter(T,Q, 14, beta(:,ii)) 
ttt=xlabel('T'); ttt.FontSize=12;
ylabel('\tau')
title(sprintf('\\beta^o_%d',ii),'FontSize', 12)
uuu=colorbar;
color_limits=[color_limits;uuu.Limits]; uuu.TickLength=0.05;
end

colormap default
for ii = 1:9
subplot(3,3,ii)
scatter(T,Q, 14, theta2(:,ii)) 
ttt=xlabel('T'); ttt.FontSize=12;
ylabel('\tau')
title(sprintf('\\beta_%d',ii),'FontSize', 12)
caxis(color_limits(ii,:))
colorbar
end


%% Run the knn admm algorithm determing  lambda using BIC

lam_set = 0.5:0.1:1.2;

theta2_set = cell(1,length(lam_set)); iii=0; lik_set =[]; sparsity_set=[];
for lambda = lam_set;
    iii=iii+1;
[theta2, theta_set2] = VC_qt_knn_admm(X,Q,T,Y,5,lambda,500,0.1,0.5);
theta2_set{iii}=theta2;
[values_lik] = likelihood_knn(theta2,X,Q,T,Y,5,lambda,500,0.1,0.5); 
lik_set=[lik_set,values_lik];

theta2_round = round(theta2,2);
sparsity=[];
for ii=1:9
    sparsity=[sparsity,length(unique(theta2_round(:,ii)))];
end
sparsity_set=[sparsity_set;sparsity];
end

%% BIC criterion
BIC_criterion_values = lik_set' + sum(sparsity_set,2) *log(n)/n;
min_BIC = min(BIC_criterion_values);
min_ind = find(BIC_criterion_values==min_BIC)
lambda_opt = lam_set(min_ind);
theta2_opt =  theta2_set{min_ind}; % the optimal result using BIC


colormap default
for ii = 1:9
subplot(3,3,ii)
scatter(T,Q, 14, theta2_opt(:,ii)) 
ttt=xlabel('T'); ttt.FontSize=12;
ylabel('\tau')
title(sprintf('\\beta_%d',ii),'FontSize', 12)
caxis(color_limits(ii,:))
colorbar
end



