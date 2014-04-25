function newM = interpM(M, r)
% function newM = interpM(M, r)
%
% This just simply interpolates along the columns of a matrix.  If you have a matrix where you want to do
% the rows, transpose it, then transpose the output.

newM = zeros(r*size(M,1), size(M,2));
for ii=1:size(M,2)
    newM(:,ii) = interp(M(:,ii),r);
end

    