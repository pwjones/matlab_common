function [y,q] = segments(x,L,M)
% function [y,q] = segments(x,L,M)
%
% x - data
% L - length of each segment
% M - the offset (in samples) between each segment
%
% Selection of segments for bootstrapping time series data.  
% Copied from Zoubir and Iskander, Bootstrap Techniques for Signal
% Processing, Cambridge University Press.


x = x(:);
N=length(x);
q=fix((N-L)/M)+1;
y=zeros(L,q);
for ii=1:q
    y(:,ii)=x((ii-1)*M+1:(ii-1)*M+L);
end;

