function [ah, fh] = subplot_pete(rows, cols, xlabeltxt, ylabeltxt, titletxt)

%nplot = rows*cols;
jj=1;
fh = figure;
for c=1:cols
    for r = rows:-1:1
        ah(jj) = axes('Position', [.05*c+(.8/cols*(c-1)) .9*(.05+.02*r)+((r-1)*.75/rows) .85/cols .9*(.75/rows)]);
        set(gca, 'TickDir','out');
        if r==1 xlabel(xlabeltxt); end
        if c==1 ylabel(ylabeltxt); end
        jj = jj+1;
    end
end
titleax = axes('Position', [0 0 1 1], 'Visible', 'off');
text(.1, .9, titletxt, 'FontSize', 16);
       