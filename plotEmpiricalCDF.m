function plotEmpiricalCDF(obs, step, colors, styles, varargin)
%function plotEmpiricalCDF(obs, colors, numSets)
%
% Varargin - 1) plot axis

num_sets = numel(obs);

if ~isempty(varargin)
    ah = varargin{1}; hold on;
else
    figure; ah = axes; hold on;
end
%if (length(varargin) > 1)

for ii=1:num_sets
    [f,x, flo, fup] = ecdf(obs{ii});
    if ii==1
        fi = find(f >= .5, 1,'first');
        pse = x(fi);
        plot(ah, [pse pse], [0 f(fi)], 'Linestyle', '--', 'Color', colors{ii});
        text('Parent', ah, 'Position', [pse-.1 .6], 'String', num2str(pse)); 
    end
    plot_err_poly(ah, x, f, f-flo, colors{ii}, ([1 1 1] + colors{ii})/2, .5);
%     ds = obs{ii};
%     minv = nanmin(ds); maxv = nanmax(ds);
%     valrange = floor(minv):step:ceil(maxv);
%     cnt = histc(ds, valrange);
%     ecdf = cumsum(cnt);
%     ecdf = ecdf./max(ecdf);
%     stairs(valrange, ecdf, 'Color', colors{ii}, 'LineStyle', styles{ii}, 'LineWidth', 1.5);
end

