
function slopeTrace = computeMiniSlope(miniTrace, dt, varargin)
% function slopeTrace = computeMiniSlope(miniTrace, dt, varargin)
%
% varargin is the size of the boxcar filter in number of samples
if ~isempty(varargin)
    filtWindSize = varargin{1};
else
    filtWindSize = 11; %filter sample size, done before taking derivative to reduce noise
end
    
miniTraceFilt = filtfilt(ones(1,filtWindSize)/filtWindSize,1,miniTrace);
slopeTrace = diff(miniTraceFilt)/dt;