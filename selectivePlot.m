function selectivePlot(x, Y, i, j) 
for ii = (i:j) - i + 1
    ax(ii) = subplot(j-i+1,1,ii);
    plot(x, Y(:,ii-1+i));
    grid on;
end
linkaxes(ax, 'xy');