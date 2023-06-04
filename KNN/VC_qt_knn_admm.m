
function [theta, theta_set] = VC_qt_knn_admm(X,Q,T,y,K,lambda,Niter,tolerance,eta)
% VC_qt_knn_admm(X,Q,T,y,5,lambda,50,0.5,1)
% X: n by p;    Q: 1 by n;   y: n;  T:1 by n  
n=size(X,1); p=size(X,2);

% initial Z
Z0 = zeros(n,p);
Z=Z0;

% initial u
u = zeros(size(Z0));

% initial theta
theta = Z0;

error = 100;
itr = 0;


theta_set = cell(1);
while error>tolerance & itr <= Niter
    itr = itr+1;
    Z0 = Z;
    theta0 = theta;

% update theta
parfor i=1:n
    xi=X(i,:); zi=Z(i,:); ui=u(i,:); 
    ci = sum(xi .* (zi-ui)) - y(i);
    if ci > (1-Q(i))* sum(xi.*xi) / eta
        thetai = zi - ui + (Q(i)-1)*xi/eta; cases=1
    end
    
    if ci < -Q(i) * sum(xi.*xi) / eta
        thetai = zi-ui + Q(i)*xi/eta; cases=2
    end
    
    if ci <= (1-Q(i))* sum(xi.*xi) / eta & ci >= -Q(i) * sum(xi.*xi) / eta
        vi = -eta * (sum(xi.*(zi-ui))-y(i))/sum(xi.*xi) ;
        thetai = zi-ui+vi*xi/eta; cases=3
    end
    theta(i,:) = thetai;    
end

% update z
Xj = [Q; T];   % quantile and index
parfor j=1:p
    yj = theta(:,j) + u(:,j);
    zj = knnfl(Xj, yj', K, lambda/eta);
    Z(:,j) = zj;
end

% update u
u = u + eta*(theta - Z);

error=max(max(abs(theta - theta0)));
[error, sqrt(sum(sum((theta - theta0).^2)))]
theta_set{itr} = theta;

end
end
