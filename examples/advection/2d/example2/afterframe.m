s = 1e-2;
axis([-0.5-s 0.5+s -0.5-s 0.5+s])
daspect([1 1 1]);

% set(pout,'edgecolor','none');

%colormap(parula);
yrbcolormap;
caxis([0,1]);

fprintf('qmin = %24.16f\n',qmin);
fprintf('qmax = %24.16f\n',qmax);


view(2);

NoQuery = 0;
prt = true;
if (prt)
    set(pout,'edgecolor','none');
    dpi = 2^7;
    figsize = [4,4];
    prefix = 'plot_nomesh';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
end

shg

clear afterframe
