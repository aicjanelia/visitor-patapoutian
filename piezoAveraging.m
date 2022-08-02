%% INSTRUCTIONS
% (1) First follow instructions to run piezoSegment_beadsRemoved
%       Note: want the _particles.mat file output by this script
% (2) Fill in the user-parameters below and run the script
% (3) The script outputs several figures for inspection. It additionally
%       generates txt files that can be loaded into PeakSelector to render
%       the super particle with brightness, etc. taken into consideration.

clc, clear, close all
%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Filenames & Directories
dataDir = 'C:\Users\RMLee\Dropbox (HHMI)\Patapoutian_iPALM_data_analyzed\iPALM_data_processed_EM\21.11.02-5\Run1-561\';
saveTag = 'Run1-561_c123_sum_X14_processed_overlay_Fiducial_transform_complete_IDL_unZ_175-285_ASCII_beadsRemoved_rRemoveX3_rRemoveY7';

% dataDir & saveTag will be used in the file name used when saving figures, etc.
spaPATH = 'C:\Users\RMLee\Documents\GitHub\smlm_datafusion3d\'; % Where is the Heydarian, et al. repository?
spaBUILD = 'C:\Users\RMLee\Documents\GitHub\Fusion3DBuild\'; % Where was the Heydarian, et al. repository built?

%%%%% Imaging Parameters
xyScale = 133.33; % nm/pix in the txt file (which is different than the rendered image!)
filtZ = 0; % If filtZ, will remove
zThresh = 75; % nm, only used if filtZ

%%%%% Scale Sweep
Nsweep = 25;

%%%%% Heydarian, et al Step 1
scale = 5; % Scaling used in the all-to-all registration, see Step 1 for more comments

%%%%% Heydarian, et al Step 2
nIterations = 5; % number of lie-algebra avg and consistency check % iterations
flagVisualizeSijHist = 1; % show S_ij histogram (boolean) 
threshold = 1; % consistency check threshold.

%%%%% Heydarian, et al Step 3
iter = 5; % number of bootstrapping iterations


%% Set up paths to dependencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath % The SPA code is sensitive to the order that folders appear in the path :|

% SPA code from Heydarian, et al.
if ~strcmp(spaPATH(end),filesep) % Make sure directory is formatted properly
    dataDir = [spaPATH filesep];
end
addpath(genpath([spaPATH 'figtree\']));
addpath(genpath([spaPATH 'expdist\']));
addpath(genpath([spaPATH 'gausstransform\']));
addpath(genpath([spaPATH 'MATLAB\']));
addpath(genpath([spaPATH 'test\']));

if ~strcmp(spaBUILD(end),filesep) % Make sure directory is formatted properly
    dataDir = [spaBUILD filesep];
end
addpath(genpath([spaBUILD 'Debug\mex']));
addpath(genpath([spaBUILD 'figtree\Debug']));

% CPU/GPU settings (CPU = 0, GPU = 1)
USE_GPU_GAUSSTRANSFORM = 1;
USE_GPU_EXPDIST = 1;

clear RESTOREDEFAULTPATH_EXECUTED

%% Load the data
if ~strcmp(dataDir(end),filesep) % Make sure directory is formatted properly
    dataDir = [dataDir filesep];
end

load([dataDir saveTag '_particles.mat'])

% Optionally filtering
if filtZ

    for kk = 1:length(particles)

        xyz = particles{kk}.coords;
        z = xyz(:,3);
        zm = mean(z(:));
        badrows = abs(z-zm) > zThresh;
%         disp([sum(badrows) size(xyz,1)])
        disp([num2str(sum(badrows)) ' of ' num2str(size(xyz,1)) ' removed'])
        particles{kk}.coords = xyz(~badrows,:);

        xyz = subParticles{kk}.points;
        z = xyz(:,3);
        zm = mean(z(:));
        badrows = abs(z-zm) > zThresh;
        subParticles{kk}.points = xyz(~badrows,:);
        subParticles{kk}.sigma = subParticles{kk}.sigma(~badrows,:);
        subParticles{kk}.txtData = subParticles{kk}.txtData(~badrows,:);
        subParticles{kk}.image = subParticles{kk}.image(~badrows,:);

    end

    saveTag = [saveTag '_filtZ' num2str(zThresh) ];

end

%% Scale sweep (Heydarian, et al)
tic

[optimal_scale,scales_vec,cost_log,idxP] = scale_sweep(particles,Nsweep,1);
saveas(gcf,[dataDir saveTag '_ScaleSweep_50particles_linspace0001-50-30.png'])
save([dataDir saveTag '_ScaleSweep_50particles_linspace0001-50-30.mat'],'optimal_scale','scales_vec','cost_log','idxP','particles')

tScaleSweep = toc;

%% Heydarian, et al Step 1
% all-to-all registration
% tic
% when the particles have a prefered orientation, like NPCs that lie in the
% cell membrane, it is recommanded to do in-plane initialization to save 
% computational time for example initAng = 1 (see pairFitting3D.m line 68). 
initAng = 'grid_72.qua';

% For NPC particle scale = 0.1 (10 nm) is the optimal choice. In other
% cases, it is recommanded to use the scale_sweep() function to find the
% optimal value. This needs to be done once for a structure (see Online
% Methods and demo2.m)

disp('all2all registration started!');
disp(datestr(now))
[RR, I] = all2all3D(subParticles, scale, initAng, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);
disp(datestr(now))
disp('all2all registration finished!');

%% Heydarian, et al Step 2
% iterations of:
% 2-1 lie-algebra averaging of relative transformation
% 2-2 consistency check
% 2-3 constructing the data-driven template

Nparticles = length(subParticles);
[initAlignedParticles, sup] = relative2absolute(subParticles, RR, I, Nparticles, ...
                                                nIterations, threshold, 1);

visualizeSMLM3D(sup,scale/2, 1); % scale/2 seems to match demo script, but they didn't list it that way...
xlabel('x'),ylabel('y'),zlabel('z')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'FontSize',16)
set(gcf, 'InvertHardcopy', 'off')
saveas(gcf,[dataDir 'initAlignedParticles_' saveTag '.png'],'png')
saveas(gcf,[dataDir 'initAlignedParticles_' saveTag '.fig'],'fig')

%% Heydarian, et al Step 3
% bootstrapping with imposing symmetry prior knowledge
USE_SYMMETRY = 3;   % flag for imposing symmetry prio knowledge (n-fold = #) <-- THIS DOES NOT WORK YET
M1 = [];            % not implemented
[superParticleWithPK, ~] = one2all3D(initAlignedParticles, iter, M1, [dataDir 'withSym_' saveTag '_'], sup, USE_SYMMETRY, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);

% % STEP 3
% bootstrapping withOUT imposing symmetry prior knowledge
USE_SYMMETRY = 0;   % flag for imposing symmetry prio knowledge
M1=[];              % not implemented
[superParticleWithoutPK, ~] = one2all3D(initAlignedParticles, iter, M1, [dataDir 'withoutSym_' saveTag '_'], sup, USE_SYMMETRY, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);

% tSteps123 = toc;

%% Convert to peakselector txt option

toTxt = [];
for kk = 1:length(initAlignedParticles)

    % Newly aligned stuff
    xyz = initAlignedParticles{kk}.points;
    sigXY = initAlignedParticles{kk}.sigma(:,1);
    sigZ = initAlignedParticles{kk}.sigma(:,2);

    % Original data
    peakSelector = initAlignedParticles{kk}.txtData;
    peakSelector = peakSelector(peakSelector(:,27) == 1,:); % Only the blades

    % Merge them!
    peakSelector(:,3) = xyz(:,1)/xyScale;
    peakSelector(:,4) = xyz(:,2)/xyScale;
    peakSelector(:,35) = xyz(:,3); % put it in z
    peakSelector(:,17) = sigXY/xyScale;
    peakSelector(:,18) = sigXY/xyScale;
    peakSelector(:,36) = sigZ;

    toTxt = [toTxt ; peakSelector];

end

% Do things go better with positive coords?
toTxt(:,3) = toTxt(:,3) - min(toTxt(:,3)) + 1;
toTxt(:,4) = toTxt(:,4) - min(toTxt(:,4)) + 1;
toTxt(:,35) = toTxt(:,35) - min(toTxt(:,35)) + 1;

writematrix(toTxt,[dataDir 'PeakSelectorFormat_' saveTag '_initAlignedParticles_scale' num2str(scale) '.txt'],'delimiter','\t')


%% Convert to peakselector txt option - bootstrapped!

x = superParticleWithoutPK{1,6}(:,1)/xyScale;
y = superParticleWithoutPK{1,6}(:,2)/xyScale;
z = superParticleWithoutPK{1,6}(:,3);

toTxt(:,3) = x;
toTxt(:,4) = y;
toTxt(:,35) = z;

%Do things go better with positive coords?
toTxt(:,3) = toTxt(:,3) - min(toTxt(:,3)) + 1;
toTxt(:,4) = toTxt(:,4) - min(toTxt(:,4)) + 1;
toTxt(:,35) = toTxt(:,35) - min(toTxt(:,35)) + 1;

writematrix(toTxt,[dataDir 'PeakSelectorFormat_' saveTag '_bootstrappedParticles_' num2str(iter) 'iterBoot_scale' num2str(scale) '.txt'],'delimiter','\t')

%% Convert to peakselector txt option - bootstrapped with symmetry

x = superParticleWithPK{1,6}(:,1)/xyScale;
y = superParticleWithPK{1,6}(:,2)/xyScale;
z = superParticleWithPK{1,6}(:,3);

toTxt(:,3) = x;
toTxt(:,4) = y;
toTxt(:,35) = z;

%Do things go better with positive coords?
toTxt(:,3) = toTxt(:,3) - min(toTxt(:,3)) + 1;
toTxt(:,4) = toTxt(:,4) - min(toTxt(:,4)) + 1;
toTxt(:,35) = toTxt(:,35) - min(toTxt(:,35)) + 1;

writematrix(toTxt,[dataDir 'PeakSelectorFormat_' saveTag '_bootstrappedSymmetryParticles_' num2str(iter) 'iterBoot_scale' num2str(scale) '.txt'],'delimiter','\t')

%% Visualize the Results - Without Symmetry

visualizeSMLM3D(superParticleWithoutPK{1,6},scale/2, 1);
xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'FontSize',16)
set(gcf, 'InvertHardcopy', 'off')
saveas(gcf,[dataDir 'superParticleNoSym_' saveTag '.png'],'png')
saveas(gcf,[dataDir 'superParticleNoSym_' saveTag '.fig'],'fig')

%% Visualize the Results - With Symmetry

visualizeSMLM3D(superParticleWithPK{1,6},scale/2, 1);
xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'FontSize',16)
set(gcf, 'InvertHardcopy', 'off')
saveas(gcf,[dataDir 'superParticleWithSym_' saveTag '.png'],'png')
saveas(gcf,[dataDir 'superParticleWithSym_' saveTag '.fig'],'fig')

%% Plot some example aligned particles
figure;
N = 10;
for kk = 1:N
    plot3(initAlignedParticles{kk}.points(:,1)+100*kk,initAlignedParticles{kk}.points(:,2)+100*kk,initAlignedParticles{kk}.points(:,3)+100*kk,'.')
    hold on
end
for kk = N+1:2*N
    plot3(initAlignedParticles{kk}.points(:,1)+100*kk+50,initAlignedParticles{kk}.points(:,2)+100*(kk-5),initAlignedParticles{kk}.points(:,3)+100*kk,'.')
end
hold off

% xlabel('X','FontSize',20)
% ylabel('Y','FontSize',20)
% zlabel('Z','FontSize',20)

set(gca,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{});
set(gca,'DataAspectRatio',[1 1 1])
view(2)
axis off
saveas(gcf,[dataDir 'AlignedParticles_' saveTag '.png'],'png')
saveas(gcf,[dataDir 'AlignedParticles_' saveTag '.fig'],'fig')

%% Figure: Localizations per particle
L = NaN*ones(length(subParticles),1);
for kk = 1:length(subParticles)
    L(kk) = size(subParticles{kk}.points,1);
end

figure(108)
histogram(L,0:5:90)
xlabel('Number of Localization Points','FontSize',20)
set(gca,'FontSize',20)
ylabel('Particle Count','FontSize',20)
box off
saveas(gcf,[dataDir saveTag '_HistogramParticleLocalizationNumber.png'],'png')

%% Save Everything
save([dataDir saveTag '_piezoAveragingWorkspace.mat'],'-v7.3')
