function [ x ] = graphtv( y, edges1,edges2, lambda,varargin )
%         z = graphtv((theta + u)', edgeunique(:,1), edgeunique(:,2), lambda/R)';
%GRAPHTV Summary of this function goes here
%   Detailed explanation goes here

if nargin==5
    weights=varargin{1};
    x = graphtv_mex(y,int32(edges1),int32(edges2),lambda,weights);
else
    x=graphtv_mex(y,int32(edges1),int32(edges2),lambda);
end

end

