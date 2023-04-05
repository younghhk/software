function [values_lik] = likelihood_knn(ini_beta,X,Q,T,y,K,lambda,Niter,tolerance,eta)
% [values_lik] =         likelihood_knn(theta2,X,Q,T,YY,KKK,lambda,100,0.1,1); 

       % [theta, theta_set VC_qt_knn_admm(X,Q,T,Y,5,lambda,150,0.2,0.5);

% X: n by p;    Q: 1 by n;   y: 1 by n;  T:1 by n  

    n = length(y); p = size(X,2);
    theta = y;
    z = y;
    u = zeros(1,n);
    QT=[Q;T];  %needed for distance measure
    y=y'; Q=Q'; T=T';
    
    Id1 = nearestneighbour(QT, 'num', K)';
    Id2 = repmat(1:n, 1, K);
    edge = [Id1(:) Id2(:)];
    data = repmat([0 0], [length(edge),1]) - edge;
    dist = sqrt(data(:,1).^2 + data(:,2).^2);
    [c,ia,ib] = unique(dist);
    edgeunique = edge(ia,:);
    edgeunique = sort(edgeunique, 2);
    
    %%
    G = graph(edgeunique(:,1), edgeunique(:,2), ones(1,size(edgeunique,1)));
    Incidence = incidence(G);
   % Laplacian = laplacian(G);
    bins = conncomp(G);
    
    row_order=cell(1,max(bins)); r_order=[]; row_order2 = cell(1,max(bins));
    Incidence_ordered = Incidence; n1=0;
    for ii=1:max(bins)
        row_order{ii}=find(bins==ii); r_order=[r_order,find(bins==ii)];
        row_order2{ii}=(n1+1):(n1+length(find(bins==ii)));
        n1=n1+length(find(bins==ii));
    end
    
    Incidence=full(Incidence);
    Incidence_ordered = Incidence(r_order,:);
    
    new_Incidence = [];
    for ii=1:size(Incidence_ordered,2)
        new_Incidence=[new_Incidence;find(Incidence_ordered(:,ii))'];      
    end
    new_G = graph(new_Incidence(:,1), new_Incidence(:,2), ones(1,size(new_Incidence,1)));

    
    [Tree,pred_vec] = minspantree(new_G,'Type','forest');
    Tree_Incidence = incidence(Tree);
   
    G_Edges = table2array(new_G.Edges); G_Edges=G_Edges(:,1:2);
    Tree_Edges = table2array(Tree.Edges);
    Tree_Edges=Tree_Edges(:,1:2);
    
    remaining_edges = G_Edges(find(1-ismember(G_Edges,Tree_Edges,'rows')),:);

      
    
    H_matrix = cell(1,max(bins)); H_matrix2=cell(1,max(bins)); H_matrix_real=[]; bar_H=[]; dot_H=[];
    arbitrary_loc=[];
   parfor ii=1:max(bins)
       H_matrix{ii}=Tree_Edges(find(ismember(Tree_Edges(:,1),row_order2{ii})),:);
       hmat = zeros(size(H_matrix{ii},1)+1,size(H_matrix{ii},1)+1);
       diff=min(min(H_matrix{ii}))-1;
       for jj=1:(size(hmat,1)-1)
           hmat(jj,H_matrix{ii}(jj,:)-diff)=[1,-1];        
       end  
       hmat(end,:) = 1./sqrt(length(hmat(end,:)));    invH=inv(hmat); 
       H_matrix{ii}=hmat; bar_H = blkdiag(bar_H, invH); 
       H_matrix_real = blkdiag(H_matrix_real,hmat);
       arbitrary_loc =[arbitrary_loc, size(hmat,1)];
       
       hmatii=remaining_edges(ismember(remaining_edges(:,1), row_order2{ii}),:);
       
       if size(hmatii,1)==0
           hmat2=[];
       else       
       hmat2=zeros(size(hmatii,1),size(hmat,2));
       for jj=1:(size(hmatii,1))
           hmat2(jj,hmatii(jj,:)-diff)=[1,-1];        
       end
       end
       H_matrix2{ii}=hmat2;
       
       dot_H = blkdiag(dot_H, hmat2*invH);
   end
   
   arbitrary_loc=cumsum(arbitrary_loc);
   
   arbitrary_loc2=[]; 
   for jj = 1:p
       arbitrary_loc2 = [arbitrary_loc2, arbitrary_loc + n*(jj-1)]      ;
   end
   
   
   mm=size(dot_H,1); dot_H = sparse(dot_H); 
   
   X_ordered = X(r_order,:); y_ordered = y(r_order); Q_ordered = Q(r_order); T_ordered = T(r_order);
  
   til_X = [];  % n by np
   ddot_H = []; % m by np
   H_matrix_real_full=[]; inv_H_matrix_real_full=[];
   parfor jj = 1:p
       til_X=[til_X,sparse(diag(X_ordered(:,jj))) * sparse(bar_H)];    
       ddot_H=blkdiag(ddot_H,dot_H);
       %H_matrix_real_full=blkdiag(H_matrix_real_full, inv(H_matrix_real));
       %inv_H_matrix_real_full = blkdiag(inv_H_matrix_real_full, H_matrix_real);
   end
   
   penidx = true(size(til_X,2),1);  % leave intercept unpenalized
   penidx(arbitrary_loc2) = false;
   
   %til_X=full(til_X);
       
    Z = y_ordered;
    theta = zeros(mm,1);
    beta = zeros(n*p,1);
    ini_beta=ini_beta(r_order,:);
    Gam = zeros(n,1); t_Gam =  zeros(mm*p,1);
    
theta_set = cell(1); itr=1; real_coef=0; error=10;
theta_set{1}=ini_beta;
%beta= inv_H_matrix_real_full * reshape(ini_beta,[],1);

beta=[];
parfor jj = 1:p
    beta =[beta; H_matrix_real * ini_beta(:,jj)];    
end


Z=y_ordered - til_X * beta;
theta =  ddot_H * beta;
beta_first=beta;

% Q_ordered
values_lik=quantile_values(Q_ordered, Z);






   
   