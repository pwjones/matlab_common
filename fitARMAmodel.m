function [Coef, CoefCI, e, fit] = fitARMAmodel(data, p, q)
% function [Coef, CoefCI, e, fit] = fitARMAmodel(data, p, q)
%
% A more understandable wrapper function for fitting an ARMA model to 
% a segment of data
Coef = []; CoefCI = []; e=[]; fit=[];

model = armax(data, [p 0 q 0]);