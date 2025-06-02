function Addpaths
%%% adding the paths
% Add Main
Main_path = [pwd,'/Base_Functions'];
addpath(Main_path);

% Add baseline paths
Baseline_path = [pwd,'/Baselines'];
addpath(Baseline_path);

% Add Channel models paths
Channel_path = [pwd,'/Channelmodels'];
addpath(Channel_path);

% Add Initialization paths
Channel_path = [pwd,'/Initialization'];
addpath(Channel_path);


addpath(pwd);
cd manopt;
addpath(genpath(pwd));
cd ..;

end