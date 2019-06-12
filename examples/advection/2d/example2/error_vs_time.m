function error_vs_time(varargin)

close all

if (length(varargin) == 1)
    % Timing data taken from 
    % projects/Basilisk/code/basilisk-wpa/examples/advection/2d/example2/...
    % errors.txt
    td = load('timing_data');
    e = td.e;
    t = td.t;
    labels = td.labels;
    norm2use = varargin{1};
else
    e = varargin{1};
    t = varargin{2};
    labels = varargin{3};
    norm2use = varargin{4};
end

m = length(e);

cl = {'r','b','g','m'};  % color
sy = {'.','p','s','v'};  % symbol
sz = [30,15,12,12];       % size
fc = {'none','b','g','m'};

for i = 1:m
    ph(i) = loglog(e{i}(:,norm2use+1),t{i},[cl{i},sy{i}]);
    set(ph(i),'markersize',sz(i));
    set(ph(i),'markerfacecolor',fc{i});
    hold on;
end


maxk = 0;
for k = 1:m
    maxk = max([maxk,length(e{k})]);
end

for j = 1:maxk
    ev = [];
    tv = [];
    for k = 1:m
        if (j > length(e{k}))
            el = nan;
            tl = nan;
        else
            el = e{k}(j,norm2use+1);
            tl = t{k}(j);
        end
        ev = [ev, el];
        tv = [tv, tl];
    end
    loglog(ev,tv,'k');
end

lh = legend(ph,labels);

xlabel('Error','fontsize',16);
ylabel('points.step/s','fontsize',16);
title('points.steps/time vs. error','fontsize',18);
set(gca,'fontsize',16);

if (norm2use == 1)
    xl = [9e-6, 0.1];
elseif (norm2use == 2)
    xl = [8.2e-5, 5.1e-1];
else
    xl = [1e-3, 1];
end
yl = [1e6,1e8];

axis([xl,yl])
% set(gca,'ytick',(2:2:20)*1e6);
% set(gca,'yticklabels',(2:2:16));

grid
set(gca,'gridlinestyle','-');
set(gca,'gridalpha',0.5);

shg



end