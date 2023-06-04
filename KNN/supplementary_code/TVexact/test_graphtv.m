%% Test script for general graph



img=imresize(imread('img.JPG'),0.3);
figure(1)
imshow(img);

img1=double(img)/255+0.1*randn(size(img));
figure(2)
imshow(img1);

%%

[ht,wt]=size(img1(:,:,1));
n=ht*wt;
m=(n*4-2*(ht+wt))/2;
edges1=zeros(m,1);
edges2=zeros(m,1);

D=sparse(m,n);

count=1;
for i=1:ht
    for j=1:wt
        here=(j-1)*ht+i;
        right=j*ht+i;
        down=(j-1)*ht+i+1;
        if i~=ht
            %D(count,here)=-1;
            %D(count,down)=1;
            edges1(count)=here;
            edges2(count)=down;
            count=count+1;
        end
        if j~=wt
            %D(count,here)=-1;
            %D(count,right)=1;
            edges1(count)=here;
            edges2(count)=right;
            count=count+1;
        end
    end
end

D=sparse(m,n);
D(edges1,edges2)=1;

k=3;

%% R

lambda=0.5;
%y=mu;
%k=3;
%[ x2,history ] = gtf_proj_newton1( y,D,k-1,lambda);
[ x ,history] = gtf_admm_v2( y,D,k-1,lambda,20*lambda);

figure(1)
semilogy(history','linewidth',2)
%title('k=5 L^3 prox')
xlabel('iteration')
ylabel('duality gap')
grid on;

figure(2)
hold off
plot(y,'r')
hold on;
plot(x,'b')
%ylim([0,2])
legend('Noisy', 'recovered')