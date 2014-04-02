function makeCircColorbar(start_theta, cm, vals, rev)
% function makeCircColorbar(start_theta, cm, vals)
%
% start_theta is the angle (in radians) you start plotting at, and you process
% around the circle counterclockwise.  0 is pointing to the right.

if isempty(rev)
    rev = 0;
end
c = [0 0];
rad = 10;
% Translate the colormap to a useable number of points.
nop = min(128, size(cm,1));
scale = size(cm,1)/nop;
cm_inds = floor(((0:nop-1) * scale)+1);
newcm = cm(cm_inds,:);
newvals = vals(cm_inds);

theta=linspace(0+start_theta,start_theta+(2*pi),nop+1);
theta = theta - theta(2)/2;
if rev
    theta = -1*theta;
end
rho=ones(1,nop)*rad;

nlabel = min(16,nop);
label_sub = floor(nop/nlabel);


figure; hold on;
for ii = 1:nop
    [p1(1) p1(2)] = pol2cart(theta(ii), rad);
    [p2(1) p2(2)] = pol2cart(theta(ii+1), rad);
    h = makePie(c,p1,p2,newcm(ii,:));
    [tp(1) tp(2)] = pol2cart(theta(ii) + (theta(ii+1)-theta(ii))/2, rad+rad/20);
    if mod(ii-1,label_sub)==0
        text(tp(1), tp(2), num2str(newvals(ii)));
    end
end
axis square;
colormap(newcm);


function h = makePie(c, p1, p2, fillc)

x = [c(1), p1(1), p2(1) c(1)];
y = [c(2), p1(2), p2(2) c(2)];
h = fill(x,y,fillc);
set(h, 'EdgeColor', 'none');
