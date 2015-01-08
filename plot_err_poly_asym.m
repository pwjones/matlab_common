function [lineh, patch_h] = plot_err_poly_asym (axis_h, time, y_mean, y_ci, color, patch_color, patch_alpha, varargin)
% [lineh, patch_h] = plot_err_poly (axis_h, time, v_mean, v_sd, color, patch_color, patch_alpha, varargin)
%
% For plotting a nice overlay mean with error areas around it.  Input arguments are named in a self
% explanatory way.  Varargin list is:
% 1: subsampling rate - plots every X points in the input.
% 2: Xlimits - uses the time vector to select a subset of the points to
% plot, out of range limits cause all to plot

  sub = 1;
  if (~isempty(varargin)) %the varargin is a subsampling ratio
      sub = varargin{1};
  end
  if (length(varargin) >= 2)
      tlim = varargin{2};
      inrange = find(time >= tlim(1) & time <= tlim(2));
      if (isempty(inrange)) inrange = 1:length(time); end
  else
      tlim = [time(1) time(end)];
      inrange = 1:length(time);
  end
  
  subi = intersect(1:sub:length(time), inrange);
  %need to protect the function from NaNs, otherwise it won't shade
  nonnan_v = find(~(isnan(y_mean)));
  nonnan_err = find(~(isnan(y_ci)));
  subi = intersect(intersect(subi, nonnan_v), nonnan_err); %taking the nonNaN, subsampled points
  time = time(subi);
  y_mean = y_mean(subi);
  y_ci = y_ci(subi,:);
  
  % First plot the SDs
  [x_err_poly, y_err_poly] = get_sem_poly(time, y_ci);
  
  if patch_alpha == 1
    patch_h = patch(x_err_poly,y_err_poly, patch_color,'EdgeColor', 'none', 'Parent', axis_h);
  else
    patch_h = patch(x_err_poly,y_err_poly, patch_color,'EdgeColor', 'none', 'Parent', axis_h, 'FaceAlpha', patch_alpha);
  end
  %line('Parent', gca, 'XData', x_err_poly, 'YData', y_err_poly, 'Linewidth', 2);
  %plot(x_err_poly, y_err_poly);
 
  % And finally the main line
  lineh = line('Parent', axis_h, 'XData', time, 'YData', y_mean, 'Color', color, 'Linewidth', 2);



%
% Returns your data as a polygon for bounding y at each x value with +/- y_off
%  for real nice SEM/SD plots
%
function [ret_x, ret_y] = get_sem_poly(x, y_ci)
  ret_x = zeros(2*length(x),1);
  ret_y = zeros(2*length(x),1);
%   max_off = max(y_off);
%   pad = max_off/1000;
%   y_off = y_off + pad; % This is because I think that plotting weirdness may be due to zero err
  l = length(x);

  for i=1:length(ret_x)
    if (i < l)
      ret_x(i) = x(i);
      ret_y(i) = y_ci(i,1);
    elseif (i == l || i == l+1)
      ret_x(i) = x(l);
      ret_y(i) = y_ci(l,1);
      if ( i == l + 1) ; ret_y(i) = y_ci(l,2); end
    else
      ret_x(i) = x(2*l - i + 1);
      %xi = mod(i-1,l)+1;
      %ret_x(i) = x(xi);
      ret_y(i) = y_ci(2*l - i + 1, 2);
    end
  end
 
