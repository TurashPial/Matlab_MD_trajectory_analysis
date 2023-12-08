clear
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
%%
A = cell2mat(xdata);
B = cell2mat(ydata);
t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,A(1,:),B(1,:),'-b')
ax1.XColor = 'b';
ax1.YColor = 'k';
%%
ax2 = axes(t);
plot(ax2,A(2,:),B(2,:),'-r')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.XColor = 'r';
%%
ax1.Box = 'off';
ax2.Box = 'off';