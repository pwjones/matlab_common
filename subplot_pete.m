function [ah, fh] = subplot_pete(rows, cols, xlabeltxt, ylabeltxt, titletxt)

nplot = rows*cols;
jj=1;
for r = 1:rows
    for c=1:cols
        ah(jj) = axes('Position', [.05*c+(.8/ncol*(c-1)) .9*(.05+.02*r+((r-1)*.75/nrow)) .85/ncol .9*(.75/nrow)]);
        set(gca, 'TickDir','out');
        xlabel(xlabeltxt);
        ylabel(ylabeltxt);
    end
end
titleax = axes('Position', [0 0 1 1], 'Visible', 'off');
fh = gcf;
text(.1, .9, titletxt, 'FontSize', 16);
       