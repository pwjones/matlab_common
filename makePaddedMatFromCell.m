function padMat = makePaddedMatFromCell(cell_in)
% function padMat = makePaddedMatFromCell(cell_in)
%
% Need a function that just takes a cell array containing different 
% lengths of vectors and make them into a matrix with the original
% vectors in the first dimension, and each vector along the second.
% Since the vectors are all different sizes, they are padded to the maximum
% length

d2 = length(cell_in);
lens = zeros(d2,1);
for ii = 1:d2
    lens(ii) = numel(cell_in{ii});
end
padMat = zeros(max(lens), d2);

for ii = 1:d2
    padMat(1:lens(ii), ii) = cell_in{ii};
end

