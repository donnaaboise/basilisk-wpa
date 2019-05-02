L = 4000;
b = 12.2;

plot_wpa = true;
plot_bas = false;
shift_t = 1;


if (shift_t)
    if (t < 30)
        axis([-50,50,-0.02,0.15]);
    else
        axis([-100,100,-0.02,0.08]);
    end
end

hclaw = getlegendinfo();
set(hclaw,'linewidth',1,'markersize',10);
xdata = get(hclaw,'xdata');
set(hclaw,'xdata',xdata - shift_t*(b + t));

fprintf('%20s = %g\n','max(AMRCLAW)',max(q));
fprintf('%20s = %d\n','length(AMRClaw)',length(q));
fprintf('\n');

lh_list = [hclaw];
lh_str = {'WPA (AMRCLAW)'};

if plot_wpa
    fname = sprintf('swe/t-%d',Frame);
    if exist(fname,'file')
        data_wpa = load(fname);
        t_wpa = data_wpa(1,1);
        if abs(t_wpa-t) > 1e-8
            error('Times do not match; t = %16.8f t_bas = %16.8f\n',t,t_wpa);
        end
        xwpa = data_wpa(2:end,1);
        hwpa = data_wpa(2:end,2);
        bthy = -1;
        etawpa = hwpa + bthy;
        
        hold on;
        b = 12.2;
        hwpa = plot(xwpa-shift_t*(b+t), etawpa, 'b.-','linewidth',1,'markersize',10);
        hold off;
        lh_list = [lh_list,hwpa];
        lh_str{end+1} = 'WPA (adaptive)';
        
        fprintf('%20s = %g\n','max(WPA)',max(etawpa));
        fprintf('%20s = %d\n','length(WPA)',length(etawpa));
        if (length(etawpa) == length(q))
            % fprintf('Diff (WPA) = %12.4e\n',norm(etawpa-q));
        end
        fprintf('\n');
    end
end

if plot_bas    
    fname = sprintf('../example3/green-naghdi/t-%d',Frame);
    if (exist(fname,'file'))
        data_bas = load(fname);
        xbas = data_bas(2:end,1);
        etabas = data_bas(2:end,2);
        
        hold on;
        b = 12.2;
        hbas = plot(xbas+shift_t*(b+t), etabas, 'ko-','linewidth',1,'markersize',5);
        hold off;
        lh_list = [lh_list,hbas];
        lh_str{end+1} = 'BAS (adaptive)';
        
        fprintf('%20s = %g\n','max(BAS)',max(etawpa));
        fprintf('%20s = %d\n','length(BAS)',length(etabas));
        if (length(etabas) == length(q))
            %fprintf('Diff = %12.4e\n',norm(etabas-q));
        end
    end
end    
lh = legend(lh_list,lh_str,'location','northeast');
set(lh,'fontsize',16);

% fprintf('q[1] = %24.16f\n',q(1));

title(sprintf('t = %g',t),'fontsize',16);

if Frame == 4 && shift_t
    axis([54, 64, -0.002, 0.03])
end

shg