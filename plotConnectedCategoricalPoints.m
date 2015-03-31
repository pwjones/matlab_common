function plotConnectedCategoricalPoints(ah, xdata, ydata, sig, varargin)

if isempty(ah)
    figure; ah = axes;
end

if nargin > 4
    err = varargin{1};
else
    err = NaN*zeros(size(xdata));
end

plot(ah, xdata, ydata, 'o-k', 'MarkerSize', 10); 
%hold on;

for ii = 1:length(sig(:))
    if (sig(ii))
        plot(ah, xdata(ii), ydata(ii), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
    end
end

        