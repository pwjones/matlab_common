function array = catStructArray(structarray, field, varargin)
%
% Varargin is a set of indices to subindex the values in the struct array field
if ~isempty(varargin)
    subinds = varargin{1};
    usesubinds = 1;
else
    usesubinds = 0;
end

array = [];
for ii = 1:length(structarray)
    if usesubinds
        temp = structarray(ii).(field)(subinds{ii});
    else
        temp = structarray(ii).(field);
    end
    array = cat(1, array, temp);
end