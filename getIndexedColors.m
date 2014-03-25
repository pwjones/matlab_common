function [cm, inds] = getIndexedColors(cm_name, data, varargin)
%function [cm, inds] = getIndexedColors(cm_name, data)
%
% So, this is the difficult way of getting color gradients for a
% pseudocolor plotting, but there are times when you can't use the colormap
% of the axis for what you want to plot.  This function helps get around
% that.
% varargin{1} is beta if you want to brighten the colormap in a way that
% preserves a good color granularity.

if ~isempty(varargin)
    equalize = varargin{1};
else
    equalize = 0;
end

dmin = nanmin(data);
dmax = nanmax(data);
dn = 2*length(data);
%dn = 256;
%if ~isempty(beta)
%    dn = ceil(length(data)*(1/(1-abs(beta))));
%else
%    dn = length(data);
%end
cm = colormap(eval([cm_name '(' num2str(dn) ')']));
inds = round((data - dmin)./(dmax-dmin) * (dn-1)) + 1;

if equalize % some extra steps if we want to stretch the colormap to best represent the data
    gm = colormap(eval(['gray(' num2str(dn) ')'])); 
    cm2 = histeq(inds, cm);
    ihist = hist(inds, dn)';
    icdf = cumsum(ihist);
    eq = (1:length(icdf))';
    norm_icdf = icdf*max(eq)./max(icdf);
    %figure; plot(eq, eq, eq, norm_icdf);
    cm = cm(round(norm_icdf),:);
    %scaling = round(
    %Ginds = round((size(gm2,1)-1)*gm2(:,1))+1;
    %cm = cm(Ginds,:);
    %cm = cm2;
end



