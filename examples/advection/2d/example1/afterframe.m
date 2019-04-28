s = 1e-2;
axis([-1-s 1+s -1-s 1+s])
daspect([1 1 1]);

% set(pout,'edgecolor','none');

colormap(parula);
caxis([-1,1]);

view(2);

NoQuery = 0;
prt = false;
if (prt)
  filename = framename(Frame,'adv0000','png');
  print('-dpng',filename);
end

shg

clear afterframe
