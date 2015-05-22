function [tout,rout] = rose2(varargin)
%ROSE   Angle histogram plot.
%   ROSE(THETA) plots the angle histogram for the angles in THETA.  
%   The angles in the vector THETA must be specified in radians.
%
%   ROSE(THETA,N) where N is a scalar, uses N equally spaced bins 
%   from 0 to 2*PI.  The default value for N is 20.
%
%   ROSE(THETA,X) where X is a vector, draws the histogram using the
%   bins specified in X. Where the values of x specify the center 
%   angle of each bin.
%
%   ROSE(AX,...) plots into AX instead of GCA.
%
%   H = ROSE(...) returns a vector of line handles.
%
%   [T,R] = ROSE(...) returns the vectors T and R such that 
%   POLAR(T,R) is the histogram.  No plot is drawn.
%
%   See also HIST, POLAR, COMPASS.

%   Clay M. Thompson 7-9-91
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5.14.4.4 $  $Date: 2005/04/28 19:56:53 $

[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,2,nargs,'struct'));

theta = args{1};
if nargs > 1, 
  x = args{2}; 
end

if ischar(theta)
  error(id('NonNumericInput'),'Input arguments must be numeric.');
end
theta = rem(rem(theta,2*pi)+2*pi,2*pi); % Make sure 0 <= theta <= 2*pi
if nargs==1,
  x = (0:19)*pi/10+pi/20;

elseif nargs==2,
  if ischar(x)
    error(id('NonNumericInput'),'Input arguments must be numeric.');
  end
  if length(x)==1,
    x = (0:x-1)*2*pi/x + pi/x;
  else
    x = sort(rem(x(:)',2*pi));
  end

end
if ischar(x) || ischar(theta)
  error(id('NonNumericInput'),'Input arguments must be numeric.');
end

% Determine bin edges and get histogram
edges = sort(rem([(x(2:end)+x(1:end-1))/2 (x(end)+x(1)+2*pi)/2],2*pi));
edges = [edges edges(1)+2*pi];
nn = histc(rem(theta+2*pi-edges(1),2*pi),edges-edges(1));
nn(end-1) = nn(end-1)+nn(end);
nn(end) = [];

% Form radius values for histogram triangle
if min(size(nn))==1, % Vector
  nn = nn(:); 
end
[m,n] = size(nn);
mm = 4*m;
r = zeros(mm,n);
r(2:4:mm,:) = nn;
r(3:4:mm,:) = nn;

% Form theta values for histogram triangle from triangle centers (xx)
zz = edges;

t = zeros(mm,1);
t(2:4:mm) = zz(1:m);
t(3:4:mm) = zz(2:m+1);

rlim = ceil(max(nn)/50)*50; %set the limit to the nearest 50
if nargout<2
  if ~isempty(cax)
    h = polar2(cax,t,r, rlim);
  else
    h = polar2(t,r,'', rlim);
  end

[a,b] = pol2cart(t,r);     % convert histogram line to polar coordinates
A = reshape(a,4,numel(x)); % reshape 4*N-element x vector into N columns
B = reshape(b,4,numel(x)); % reshape 4*N-element y vector into N columns
patch(A,B,[.3 .6 1])         % plot N patches based on the columns of A and B
  
  if nargout==1, tout = h; end
  return
end

if min(size(nn))==1,
  tout = t'; rout = r';
else
  tout = t; rout = r;
end

function str=id(str)
str = ['MATLAB:rose:' str];

