OutputFlag = 'ascii';
% OutputDir = './uniform';
OutputDir = '_output';

mq = 1;                       % which component(s) of q to plot
UserVariable = 1;             % set to 1 to specify a user-defined variable
UserVariableFile = 'surface_eta';       % name of m-file mapping data to q
% UserVariableFile = 'uvelocity';       % name of m-file mapping data to q
MappedGrid = 0;               % set to 1 if mapc2p.m exists for nonuniform grid
PlotStyle = setplotstyle('r.-');  % used in plot command for line color and type
MaxFrames = 1000;                % max number of frames to loop over
