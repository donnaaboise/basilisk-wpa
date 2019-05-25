function error_vs_time(e,t,labels,norm2use)

close all



% e1 = e_bas(:,1+norm2use);
% t1 = t_bas;
% 
% e2 = e_wpa(:,1+norm2use);
% t2 = t_wpa;

cl = {'r','b','g'};  % color
sy = {'.','p','s'};  % symbol
sz = [30,15,12];       % size
fc = {'none','b','g'};

for i = 1:length(e)
    ph(i) = loglog(e{i}(:,norm2use+1),t{i},[cl{i},sy{i}]);
    set(ph(i),'markersize',sz(i));
    set(ph(i),'markerfacecolor',fc{i});
    hold on;
end


maxk = 100;
for i = 1:3
    if (length(e{i}) < maxk)
        maxk = length(e{i});
    end
end

for k = 1:maxk
    ev = [e{1}(k,norm2use+1),e{2}(k,norm2use+1),e{3}(k,norm2use+1)];
    tv = [t{1}(k),t{2}(k),t{3}(k)];
    loglog(ev,tv,'k');
end

lh = legend(ph,labels);

xlabel('Error','fontsize',16);
ylabel('points.step/s','fontsize',16);
title('points.steps/time vs. error','fontsize',18);
set(gca,'fontsize',16);

xl = [1e-4, 0.1];
yl = [1e6,1e8];

axis([xl,yl])
% set(gca,'ytick',(2:2:20)*1e6);
% set(gca,'yticklabels',(2:2:16));

grid
set(gca,'gridlinestyle','-');
set(gca,'gridalpha',0.5);


shg



end