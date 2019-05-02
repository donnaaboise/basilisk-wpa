s = 1e-2;
axis([-0.5-s 0.5+s -0.5-s 0.5+s])
daspect([1 1 1]);

% set(pout,'edgecolor','none');

fprintf('qmin = %20.16f\n',qmin);
fprintf('qmax = %20.16f\n',qmax);

colormap(parula);
if (Frame == 25)
    % Use qmin, qmax from ForestClaw
    caxis([0.0702192686928361, 0.1691585563679505]);
else
    caxis([qmin,qmax]);
end
ss = 1e-12;
caxis([0.1-ss, 0.1+ss]);

NoQuery = 0;
prt = true;
if (prt)
    set(pout,'edgecolor','k');
    maxlevel = 10;
    dpi = 2^7;    % fix at 128
    figsize = (2^maxlevel/dpi)*[1,1];
    prefix = 'plot_wpa';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
end

shg 


%clear afterframe
