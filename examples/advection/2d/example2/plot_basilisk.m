function plot_basilisk(dir,frames)

maxframes = 1000;

if (nargin < 2)
    frames = 1:maxframes;
end

done = false;
frameno = 0;
while(~done)
    t = frames(frameno+1);
    
    plot_data(dir,frameno);
    
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

function plot_data(dir,frameno)

fname = sprintf('%s/fort.t%04d',dir,frameno);
if ~exist(fname,'file')
    fprintf('\n');
    fprintf('****** File %s does not exist\n',fname);
    fprintf('\n');
    pout = 0;
    return;
end
fp = fopen(fname,'r');
data = fscanf(fp,'%g',1); fscanf(fp,'%s',1);
meqn = fscanf(fp,'%g',1); 
fclose(fp);
t = data(1,1);

fprintf('\n');
fprintf('Frame %d at t = %g\n',frameno, t);
fprintf('Reading data from %s\n\n',fname);

fname = sprintf('%s/fort.f%04d',dir,frameno);
fpf = fopen(fname,'r');
F = fread(fpf,[5,Inf],'int')';
fclose(fpf);

fname = sprintf('%s/fort.c%04d',dir,frameno);
fpc = fopen(fname,'r');
C = fread(fpc,[meqn,Inf],'double')';
fclose(fpc);

fname = sprintf('%s/fort.v%04d',dir,frameno);
fpv = fopen(fname,'r');
V = fread(fpv,[2,Inf],'double')';
fclose(fpv);

pout = patch('Faces',F,'Vertices',V);
set(pout,'CData',C);
set(pout,'Facecolor','flat');

title(sprintf('t = %g',t),'fontsize',16);

Frame = frameno;
qmin = min(C(:));
qmax = max(C(:));
assignin('caller','qmin',qmin);
assignin('caller','qmax',qmax);
assignin('caller','pout',pout);
assignin('caller','Frame',Frame);

if exist('afterframe.m','file')
    afterframe;
end

end