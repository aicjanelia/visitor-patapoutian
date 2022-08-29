% In comparison to batch_piezoAveraging, this script creates an
% intermediate file during all2all3D which allows progress to be
% interupted and restarted more easily.
%%%%%% INSTRUCTIONS
% (1) First follow instructions to run piezoSegment_beadsRemoved
%       Note: want the _particles.mat file output by this script
% (2) Fill in the user parameters below and run the script. The variable
%       directories allows this script to loop through multiple experiments,
%       analyzing each one in turn.
% (3) The main function piezoAveraging outputs several figures for 
%       inspection. It additionally generates txt files that can be loaded 
%       into PeakSelector to render the super particle with brightness,
%       etc. taken into consideration.

clc, clear, close all

%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List all directories you would like to process
directories = {'fullpath\todataset1\',...
    'fullpath\todataset2\',...
    'fullpath\todataset3\'};

% This script is currently set to sequentially run BOTH z-filtering and no filtering
zThresh = 75; % nm, only used if filtZ

% Directories to the dependencies
spaPATH = 'fullpath\to_smlm_datafusion3d_iPALM_repository\'; % Where is the smlm_datafusion3d_iPALM repository?
spaBUILD = [spaPATH 'build\']; % Where was the repository built?

% Overwite exisiting analysis files?
overwrite = 0; % Set to true to repeat already completed analysis


%% Set up paths to dependencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath % The SPA code is sensitive to the order that folders appear in the path :|

% SPA code from Heydarian, et al.
if ~strcmp(spaPATH(end),filesep) % Make sure directory is formatted properly
    spaPATH = [spaPATH filesep];
end
addpath(genpath([spaPATH 'figtree\']));
addpath(genpath([spaPATH 'expdist\']));
addpath(genpath([spaPATH 'gausstransform\']));
addpath(genpath([spaPATH 'MATLAB\']));
addpath(genpath([spaPATH 'test\']));

if ~strcmp(spaBUILD(end),filesep) % Make sure directory is formatted properly
    spaPATH = [spaBUILD filesep];
end
addpath(genpath([spaBUILD 'Release\mex']));
addpath(genpath([spaBUILD 'figtree\Release']));

%% LOOP THROUGH DIRECTORIES TO RUN AVERAGING

for kk = 1:length(directories)

    dataDir = directories{kk};
    disp(['Processing directory ' num2str(kk) '/' num2str(length(directories)) ': ' dataDir])

    saveTag = dir([dataDir '*_particles.mat']);
    if length(saveTag) > 1
        warning('Using first found *particles.mat file')
    elseif isempty(saveTag)
        warning(['No *particles.mat file found. ' 10 dataDir ' will be skipped'])
        continue
    end
    saveTag = saveTag(1).name(1:end-14);

     % Run with z-filtering
    filtZ = 1;
    if ~exist([dataDir saveTag '_filtZ' num2str(zThresh) '_piezoAveragingWorkspace.mat'],'file') || overwrite
        piezoAveragingIntermediate(dataDir,saveTag,filtZ,zThresh,0)
        disp(' ')
        disp(' ')
    else
        disp(['Z-filtering already exisits for ' dataDir])
        disp([saveTag '_filtZ' num2str(zThresh) '_piezoAveragingWorkspace.mat already exists.'])
        disp(' ')
        disp(' ')
    end

    % Run without z-filtering
    filtZ = 0;
    if ~exist([dataDir saveTag '_piezoAveragingWorkspace.mat'],'file') || overwrite
        piezoAveragingIntermediate(dataDir,saveTag,filtZ,zThresh,0)
        disp(' ')
        disp(' ')
    else
        disp(['Unfiltered analysis already exisits for ' dataDir])
        disp([saveTag '_piezoAveragingWorkspace.mat already exists.'])
        disp(' ')
        disp(' ')
    end

end

%% WRAP UP
disp(' ')
disp('BATCH SCRIPT COMPLETE')