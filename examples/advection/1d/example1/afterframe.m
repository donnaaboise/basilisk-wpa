axis([-1,1,-0.2,1.5]);

plot_exact = true;
plot_wpa = true;

fprintf('%20s %g\n','max(q)',max(q));
fprintf('%20s %d\n','length(uniform)',length(q));
fprintf('\n');

hclaw = getlegendinfo();
set(hclaw,'markersize',5);
lh_list = [hclaw];
lh_str = {'WPA (uniform)'};

if plot_exact
    qinit = @(x) sin(2*pi*x);
    qinit = @(x) abs(x) < 0.25;
    ubar = 0.5;

    xf = linspace(-1,1,bitshift(1,16));
    xt = mod(xf+1-ubar*t,2) - 1;
    qf = qinit(xt);

    hold on;
    hf = plot(xf,qf, 'k-','linewidth',1);
    hold off;
    lh_list = [lh_list,hf];
    lh_str{end+1} = 'Exact';
end



if plot_wpa
    fname = sprintf('adv/t-%d',Frame);
    if exist(fname,'file')
        data_wpa = load(fname);
        t_wpa = data_wpa(1,1);
        if abs(t_wpa-t) > 1e-8
            error('Times do not match; t = %16.8f t_bas = %16.8f\n',t,t_wpa);
        end
        xwpa = data_wpa(2:end,1);
        qwpa = data_wpa(2:end,2);
        
        hold on;
        hwpa = plot(xwpa, qwpa, 'b.-','linewidth',1,'markersize',10);
        hold off;
        lh_list = [lh_list,hwpa];
        lh_str{end+1} = 'WPA (adaptive)';
        
        fprintf('%20s = %g\n','max(WPA)',max(qwpa));
        fprintf('%20s = %d\n','length(WPA)',length(qwpa));
        fprintf('\n');
    end
end

lh = legend(lh_list,lh_str,'location','southeast');
set(lh,'fontsize',16);

title(sprintf('t = %g',t),'fontsize',16);

shg