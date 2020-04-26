function [ah, fh] = subplot_pete(rows, cols, xlabeltxt, ylabeltxt, titletxt)

%nplot = rows*cols;
jj=1;
fh = figure;
rowh = .8/rows;
colw = .9/cols;
for c=1:cols
    for r = rows:-1:1
        ah(jj) = axes('Position', [(colw*(c-1)+.09) .09+rowh*(r-1) .8*colw .78*rowh]);
        set(gca, 'TickDir','out');
        if r==1 xlabel(xlabeltxt); end
        if c==1 ylabel(ylabeltxt); end
        jj = jj+1;
    end
end
titleax = axes('Position', [0 0 1 1], 'Visible', 'off');
text(.2, .95, titletxt, 'FontSize', 16);
       