function plot_errors(e,idx,compare)

if (nargin < 3)
    compare = false;
end

lh = zeros(3,1);
lstr_basic = {'1-norm','2-norm','inf-norm'};
lstr = lstr_basic;
if (~compare)
    close all
    ladd = 0;
else
    lh_compare = findobj('Tag','compare');
    if (~isempty(lh_compare))
        lh = lh_compare;
        o = findobj('Tag','legend');
        lstr = o.String;
        ladd = 3;
    else
        ladd = 0;
    end
end



Nvec = e(:,1);
n = length(Nvec);

c = {'r','b','g'};
if (nargin == 1)
    idx = 1:n;
end
for i = 1:3
    p = polyfit(log(Nvec(idx)),log(e(idx,i+1)),1);
    lh(i+ladd) = loglog(Nvec,exp(polyval(p,log(Nvec))),c{i},'linewidth',2);     
    lstr{i+ladd} = sprintf('%10s (rate = %6.4f)',lstr_basic{i},-p(1));
    hold on;
end
mh = loglog(Nvec,e(:,2:4),'k.','markersize',17);
hold on;

if (~compare)
    i1 = idx(1);
    i2 = idx(end);
    yl = ylim;
    plot([Nvec(i1) Nvec(i1)],yl,'k--'); 
    plot([Nvec(i2) Nvec(i2)],yl,'k--');
    set(lh,'Tag','compare');
else
    n = length(lh);
    set(lh(end-2:end),'linestyle','--');
    set(mh,'marker','p','markersize',12);
end

legend(lh,lstr,'location','southwest');

xlabel('N','fontsize',16);
ylabel('Error','fontsize',16);
title('Error','fontsize',18);

set(gca,'xtick',Nvec);

p0 = log2(Nvec(1));
p1 = log2(Nvec(end));
xlim([2^(p0-0.5), 2^(p1+0.5)]);

set(gca,'fontsize',16);

if (size(e,2) == 5)
    idx = [2:5];
else
    idx = 2:4;
end
conv_rates = log(e(1:end-1,idx)./e(2:end,idx))/log(2);
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
    


end