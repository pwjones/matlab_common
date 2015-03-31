function eHandles = addErrorBarsAsym(axisHandle, x, y, lb, ub, color, varargin)
% function eHandles = addErrorBars(axisHandle, x, y, e, color, varargin)
%
% This function just adds error bars to a graph, like the built-in
% MATLAB function except this returns the handles of the error bars
% and is fully compatible with code that uses primitive plotting 
% tools because it doesn't do anything else.
%
% varargin accepts a length of the hat added to the error bars
x = x(:);
y = y(:);
lb = lb(:);
ub = ub(:);

hatlen = 1;
if (~isempty(varargin))
     hatlen = varargin{1};
end
eHandles = zeros(length(x)*3,1);
for i = 1:length(x)
    handlei = 3*(i-1)+1;
    top = ub(i);
    bottom =  lb(i);
    ydata = [top bottom];
    xdata = [x(i) x(i)];
    eHandles(handlei) = line('Parent', axisHandle, 'XData', xdata, 'YData', ydata, 'Color', color, 'Tag', 'error');
    xhat = [x(i)-hatlen x(i)+hatlen];
    ytophat = [top top];
    ybottomhat = [bottom bottom];
    eHandles(handlei+1) = line('Parent', axisHandle, 'XData', xhat, 'YData', ytophat, 'Color', color, 'Tag', 'error');
    eHandles(handlei+2) = line('Parent', axisHandle, 'XData', xhat, 'YData', ybottomhat, 'Color', color, 'Tag', 'error');
    
end
    
    
    