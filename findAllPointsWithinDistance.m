function pts = findAllPointsWithinDistance(path, R)

% First idea: for each square of side length 2R, the points will be
% enclosed

nCircSamp = R*8;
theta = linspace(0,2*pi, nCircSamp+1);
pts = [];
for jj = 1:size(path, 1)
    for ii = 1:length(theta)
        [x,y] = pol2cart(theta(ii), 1:R);
        pts_temp = [x', y'] + repmat(path(jj,:), length(x), 1);
        pts = cat(1,pts, round(pts_temp));
        pts = unique(pts, 'rows');
    end
end

%figure;
%hold on;
%plot(pts(:,1), pts(:,2), 'r.');
%plot(path(:,1), path(:,2), 'ko', 'MarkerFaceColor', 'k');


