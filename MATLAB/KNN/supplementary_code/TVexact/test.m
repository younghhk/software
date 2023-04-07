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

%%



img2=img1;
%  img2(:,:,1)=tvexact(img1(:,:,1),0.5);
%  img2(:,:,2)=tvexact(img1(:,:,2),0.5);
%  img2(:,:,3)=tvexact(img1(:,:,3),0.5);

tic
 lambda=0.1;
 for i=1:3
 y=img1(:,:,i);
 y=y(:);
 [ x ] = graphtv( y, edges1,edges2, lambda);
img2(:,:,i)=reshape(x,[ht,wt]);
 end
timeelapsed=toc
%img=rgb2gray(img);

figure(3)
imshow(img2)

