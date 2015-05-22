function M = cell2mat_2Dunequal(input)
% function cell2mat_2Dunequal(cell)
%
% Utility function for the case of having a cell array of 1D vectors of
% uneven lengths. Creates an output matrix padded by NaNs of 
% size (max(numel(input{})), numel(input))
lens = zeros(numel(input),1);
for ii = 1:numel(input)
    lens(ii) = numel(input{ii});
end
M = NaN*zeros(max(lens), numel(input));
for ii = 1:numel(input)
    M(1:numel(input{ii}), ii) = input{ii}(:);
end

