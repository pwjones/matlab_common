function fieldNum = findh5FieldNumber(infoStruct, Name)
%
%
% Just a for loop that goes through and looks for a field name that matches
% the string

N = length(infoStruct);
fieldNum = [];
for ii = 1:N
    if(strcmp(infoStruct(ii).Name, Name));
        fieldNum = ii;
        break;
    end
end
