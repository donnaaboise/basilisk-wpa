function plot_errors(varargin)

use_saved_data = varargin{1};

if (use_saved_data)
    % Timing data taken from 
    % projects/Basilisk/code/basilisk-wpa/examples/advection/2d/example2/...
    % errors.txt
    % labels='Basilisk', 'ForestClaw', 'WPA (fwaves)', 'WPA (color)'
    run2use = varargin{2};
    idx = varargin{3};
    compare = varargin{4};
    td = load('timing_data');
    e = td.e{run2use};
    t = td.t;
    labels = td.labels;
    tstr = labels{run2use};
else    
    e = varargin{2};
    idx = varargin{3};
    compare = varargin{4};
    tstr = 'Input errors';
end

lh = zeros(3,1);
lsq = zeros(3,1);
lstr_basic = {'1-norm','2-norm','inf-norm'};
lstr = lstr_basic;
if (~compare)
    close all
    ladd = 0;
else
    lh_compare = findobj('Tag','compare');
    lsq = findobj('Tag','lsqfit');
    if (~isempty(lh_compare))        
        lh = lh_compare;
        o = findobj('Tag','legend');
        lstr = o.String;
        ladd = 3;
    else
        ladd = 0;
    end
    if (~isempty(lsq))
        delete(lsq);
    end
end



Nvec = e(:,1);
n = length(Nvec);

c = {'r','b','g','c'};
m = {'.','p','v','s'};
if (nargin == 1)
    idx = 1:n;
end
for i = 1:3
    p = polyfit(log(Nvec(idx)),log(e(idx,i+1)),1);
    if (~compare)
        lsq(i) = loglog(Nvec,exp(polyval(p,log(Nvec))),'k','linewidth',1); 
        set(lsq(i),'Tag','lsqfit');
        hold on;
    end
    lh(i+ladd) = loglog(Nvec,e(:,i+1),c{i},'linewidth',1);
    lstr{i+ladd} = sprintf('%10s (rate = %6.4f)',lstr_basic{i},-p(1));
    hold on;
end

if (~compare)
    mh = loglog(Nvec,e(:,2:4),'ko','markersize',8);
    mh_idx = loglog(Nvec(idx),e(idx,2:4),'k.','markersize',28);
else
    n = length(lh);
    set(lh(end-2:end),'linestyle','--');
    mh = loglog(Nvec,e(:,2:4),'kp','markersize',12);
    % set(mh,'markerfacecolor','r');
    mh_idx = loglog(Nvec(idx),e(idx,2:4),'kp','markersize',14);
    set(mh_idx,'markerfacecolor','k');
end    
set(lh,'Tag','compare');

legend(lh,lstr,'location','southwest');

xlabel('N','fontsize',16);
ylabel('Error','fontsize',16);
title(sprintf('Error (%s)',tstr),'fontsize',18);

set(gca,'xtick',Nvec);

p0 = log2(Nvec(1));
p1 = log2(Nvec(end));
xlim([2^(p0-0.5), 2^(p1+0.5)]);

set(gca,'fontsize',16);

if (size(e,2) == 5)
    cidx = [2:5];
else
    cidx = 2:4;
end
conv_rates = log(e(1:end-1,cidx)./e(2:end,cidx))/log(2);
fprintf('\n');
fprintf('          Convergence rates\n');
cr = [Nvec(2:end) Nvec(1:end-1) conv_rates];
if (size(e,2) == 5)
    fprintf('%s\n',double('-')*ones(1,48));
    fprintf('%5d/%5d %8.4f %8.4f %8.4f %8.4f\n',cr');
    fprintf('%s\n',double('-')*ones(1,48));
else
    fprintf('%s\n',double('-')*ones(1,39));
    fprintf('%5d/%5d %8.4f %8.4f %8.4f\n',cr');
    fprintf('%s\n',double('-')*ones(1,39));
end
    
ylim([1e-6,1e1]);


shg
end