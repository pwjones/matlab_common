function ivect = indicesFromTimes(tvect, lookup)
% function ivect = indicesFromTimes(tvect, lookup)
%
% Use a time vector to find indices of a vector of time values
% Will give the index of the first time value greater than or equal
% to each time in the lookup list.
ivect = ones(size(lookup))*-1;
for ii = 1:length(lookup)
    fi = find(tvect >= lookup(ii), 1, 'first');
    if ~isempty(fi)
        ivect(ii) = fi;
    end
end
