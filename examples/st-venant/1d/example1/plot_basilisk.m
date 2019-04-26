function plot_basilisk(dir,tvec)

maxframes = 1000;

if (nargin < 2)
    tvec = 0:maxframes;
end

done = false;
frameno = 0;
while(~done)
    t = tvec(frameno+1);
    
    success = plot_data(dir,t,frameno);
    if (success)
        afterframe;
    end
    
    str = input('Hit <return> for next plot, or type k, j or q : ','s');
    switch str
        case 'k'
            keyboard;
        case 'q'
            return;
        case 'j'
            frameno = input('Enter frame number to plot : ');  
            frameno = max([0,min([frameno,maxframes])]);
        otherwise
            frameno = frameno + 1;
    end
    
    if (frameno >= maxframes)
        done = true;
    end
    

end

end

function success = plot_data(dir,t,frameno)

fname = sprintf('%s/t-%d',dir,t);
if ~exist(fname,'file')
    fprintf('\n');
    fprintf('****** File %s does not exist\n',fname);
    fprintf('\n');
    success = false;
    return;
end
data = load(sprintf('%s/t-%d',dir,t));
t = data(1,1);

fprintf('\n');
fprintf('Frame %d at t = %g\n',frameno, t);
fprintf('Reading data from %s\n\n',fname);

x = data(2:end,1);
h = data(2:end,2);
hu = data(3:end,3);
bthy = -1;

q = h + bthy;

b = 12.2;

assignin('caller','xcenter',x);
assignin('caller','t',t);
assignin('caller','q',q);
assignin('caller','Frame',frameno);

Frame = frameno;
xcenter = x;

shift_t = 0;

hclaw = plot(x - shift_t*(b + t),q,'r.-','markersize',10,'linewidth',1);

assignin('caller','hclaw',hclaw);

title(sprintf('t = %g',t),'fontsize',16);

success = true;

end