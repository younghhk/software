%% Varying coefficient function algorithm


% Samples and indices
n=3000;
T = rand(1,n);  u =T;

% randomly choosing Q (quantile)  
Q = rand(1,n);
Q = 0.9*Q + 0.05; 

% generating data
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


%% Run the knn admm algorithm with fixed lambda and K=5
% Here we set K=5, tolerance =0.1, eta (step-size parameter) = 0.5. 
lambda = 1;
[theta2, theta_set2] = VC_qt_knn_admm(X,Q,T,Y,5,lambda, 500, 0.1, 0.5);

%% underlying and estimated coefficient plots
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



