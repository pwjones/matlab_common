function slope = fitCumTsSlope(data)

cumData = cumsum(data(:));
p = polyfit((1:length(cumData))', cumData,1);
slope = p(1);