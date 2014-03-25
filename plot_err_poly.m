function [lineh, patch_h] = plot_err_poly (axis_h, time, v_mean, v_sd, color, patch_color, patch_alpha, varargin)
% [lineh, patch_h] = plot_err_poly (axis_h, time, v_mean, v_sd, color, patch_color, patch_alpha, varargin)
%
% For plotting a nice overlay mean with error areas around it.  Input arguments are named in a self
% explanatory way.  Varargin list is:
% 1: subsampling rate - plots every X points in the input.
% 2: Xlimits - uses the time vector to select a subset of the points to
% plot, out of range limits cause all to plot

  sub = 1;
  if (length(varargin)>0 && ~isempty(varargin{1})) %the varargin is a subsampling ratio
      sub = varargin{1};
  end
  if (length(varargin) >= 2 && ~isempty(varargin{2})) % a time limit
      tlim = varargin{2};
      inrange = find(time >= tlim(1) & time <= tlim(2));
      if (isempty(inrange)) inrange = 1:length(time); end
  else
      tlim = [time(1) time(end)];
      inrange = 1:length(time);
  end
  if (length(varargin) >= 3 && ~isempty(varargin{3})) % plot option pairs
      plotOptions = varargin{3};
  else
      plotOptions = {};
  end
  
  subi = intersect(1:sub:length(time), inrange);
  %need to protect the function from NaNs, otherwise it won't shade
  nonnan_v = find(~(isnan(v_mean)));
  nonnan_err = find(~(isnan(v_sd)));
  subi = intersect(intersect(subi, nonnan_v), nonnan_err); %taking the nonNaN, subsampled points
  time = time(subi);
  v_mean = v_mean(subi);
  v_sd = v_sd(subi);
  
  % First plot the SDs
  [x_err_poly, y_err_poly] = get_sem_poly(time, v_mean, v_sd);
  
  patch_h = patch(x_err_poly,y_err_poly, patch_color,'EdgeColor', 'none', 'Parent', axis_h); % , 'FaceAlpha', patch_alpha
  %line('Parent', gca, 'XData', x_err_poly, 'YData', y_err_poly, 'Linewidth', 2);
  %plot(x_err_poly, y_err_poly);
 
  % And finally the main line
  lineh = line('Parent', axis_h, 'XData', time, 'YData', v_mean, 'Color', color, 'Linewidth', 2);
  if ~isempty(plotOptions)
    set(lineh, 'LineStyle', plotOptions{1});
  end


%
% Returns your data as a polygon for bounding y at each x value with +/- y_off
%  for real nice SEM/SD plots
%
function [ret_x, ret_y] = get_sem_poly(x, y, y_off)
  ret_x = zeros(2*length(x),1);
  ret_y = zeros(2*length(x),1);
%   max_off = max(y_off);
%   pad = max_off/1000;
%   y_off = y_off + pad; % This is because I think that plotting weirdness may be due to zero err
  l = length(x);

  for i=1:length(ret_x)
    if (i < l)
      ret_x(i) = x(i);
      ret_y(i) = y(i) - y_off(i);
    elseif (i == l || i == l+1)
      ret_x(i) = x(l);
      ret_y(i) = y(l) - y_off(l);
      if ( i == l + 1) ; ret_y(i) = y(l) + y_off(l); end
    else
      ret_x(i) = x(2*l - i + 1);
      %xi = mod(i-1,l)+1;
      %ret_x(i) = x(xi);
      ret_y(i) = y(2*l - i + 1) + y_off(2*l - i + 1);
    end
  end
 
