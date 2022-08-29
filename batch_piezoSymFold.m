%%%%%% INSTRUCTIONS
% Given an already averaged set of particles, with an exisiting
% _piezoAveragingWorkspace.mat file, this script will (1) fit the points to
% a plane, (2) rotate that plane into the x-y plane, (3) peform
% bootstrapping on the newly rotated particles, and (4) generate animations
% of each resulting super particle.

clc, clear, close all

%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List all directories you would like to process
directoriesNow = {'fullpath\todataset1\',...
    'fullpath\todataset2\',...
    'fullpath\todataset3\'};

% This script is currently set to sequentially run BOTH z-filtering and no filtering
zThreshNow = 75; % nm, only used if filtZ

% Directories to the dependencies
spaPATH = 'fullpath\to_smlm_datafusion3d_iPALM_repository\'; % Where is the smlm_datafusion3d_iPALM repository?
spaBUILD = [spaPATH 'build\']; % Where was the repository built?

% Overwite exisiting analysis files?
overwrite = 1; % Set to true to repeat already completed analysis

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

for mm = 1:length(directoriesNow)

    %%% Load the data
    dataDirNow = directoriesNow{mm};
    disp(['Processing directory ' num2str(mm) '/' num2str(length(directoriesNow)) ': ' dataDirNow])

    saveTag = dir([dataDirNow '*_particles.mat']);
    if length(saveTag) > 1
        warning('Using first found *particles.mat file')
    elseif isempty(saveTag)
        warning(['No *particles.mat file found. ' 10 dataDirNow ' will be skipped'])
        continue
    end
    saveTag = saveTag(1).name(1:end-14);

    try % This will replace saveTag, be careful!
        load([dataDirNow saveTag '_filtZ' num2str(zThreshNow) '_piezoAveragingWorkspace.mat'])
    catch err
        if mm == 1 || mm == 2 % These two were not run this way
            warning('for first two directories, loading unfiltered')
            load([dataDirNow saveTag '_piezoAveragingWorkspace.mat'])
        elseif mm > 11
            load([dataDirNow saveTag '_piezoAveragingWorkspace.mat'])
            try
                load([dataDirNow saveTag(10:end) '_piezoAveragingWorkspace.mat'])
            catch err
                disp(err)
                error('Something is wrong with the z-filt file') 
            end
        else
            disp(err)
            error('Something is wrong with the z-filt file')            
        end
    end

    if exist([dataDirNow saveTag '_piezoAveragingWorkspaceRotated.mat'],'file') && ~overwrite
        continue
    end

     %% Rotate the initially aligned particles into the xy plane
     % Based on the best fit to the super imposed particle
     % See https://stackoverflow.com/questions/9423621/3d-rotations-of-a-plane

    data = sup;
    density = mvksdensity(data,data,'Bandwidth',scale/2); %  This is the probability density
    
    figure(10)
    subplot(1,2,1)
    set(gcf,'Position',[500 275 560*2 420*2])    
    scatter3(data(:,1),data(:,2), data(:,3), density*5*10^6, density, '.')
    colormap(cool)
    xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'FontSize',16)
    
    if mm <= 11
        % Manually set the axis to a range I like
        axis(40*[-1 1 -1 1 -1 1])
        caxis([10^-5 8*10^-5])
    else
        axMin = floor(min(data(:)));
        axMax = ceil(max(data(:)));
        axis(repmat([axMin axMax],[1 3]))
    end
          
    set(gcf,'Color','white')
    title('original data')

    [fitTest,gof] = fit([sup(:,1) sup(:,2)],sup(:,3),'poly11','Weights',density);
    [xtest,ytest] = meshgrid(-40:40,-40:40);
    ztest = fitTest.p00 + fitTest.p10*xtest + fitTest.p01*ytest;
    hold on
    surf(xtest,ytest,ztest,0.6*ones(size(xtest)),'FaceAlpha',0.2)
    shading flat
    hold off
    title('sup + fitted plane')

    Nvec = [0,0,1]; % zhat
    Mvec = -[fitTest.p10,fitTest.p01,-1]; % Normal to the fitted plane

    c = dot(Mvec,Nvec)/(norm(Mvec)*norm(Nvec)); %cos theta
    axR = cross(Mvec,Nvec)/norm(cross(Mvec,Nvec)); % axis of rotation
    x = axR(1); y = axR(2); z = axR(3); % components
    s = sqrt(1-c*c); % sin theta
    C = 1-c;
    rmat = [ x*x*C+c    x*y*C-z*s  x*z*C+y*s ;...
             y*x*C+z*s  y*y*C+c    y*z*C-x*s ;...
             z*x*C-y*s  z*y*C+x*s  z*z*C+c   ]; % rotation of theta around the chose axis

    newPoints = rmat*([sup(:,1) sup(:,2) sup(:,3)]');
    newPoints = newPoints';
    fitPoints = fit([newPoints(:,1) newPoints(:,2)],newPoints(:,3),'poly11','Weights',density);

    subplot(1,2,2)
    set(gcf,'Position',[500 275 560*2 420*2])    
    scatter3(newPoints(:,1),newPoints(:,2), newPoints(:,3), density*5*10^6, density, '.')
    colormap(cool)
    xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'FontSize',16)

    zPoints = fitPoints.p00 + fitPoints.p10*xtest + fitPoints.p01*ytest;
    hold on
    surf(xtest,ytest,zPoints,0.2*ones(size(xtest)),'FaceAlpha',0.2)
    shading flat
    hold off
    title('rotated')

    if mm <= 11
        % Manually set the axis to a range I like
        axis(40*[-1 1 -1 1 -1 1])
        caxis([10^-5 8*10^-5])
    else
        axMin = floor(min(data(:)));
        axMax = ceil(max(data(:)));
        axis(repmat([axMin axMax],[1 3]))
    end

    saveas(gcf,[dataDir 'RotationPlanes_' saveTag '.png'],'png')
    saveas(gcf,[dataDir 'RotationPlanes_' saveTag '.fig'],'fig')

    %%%% Now apply the rotation to each of the intially aligned particles
    initAlignedParticlesRot = initAlignedParticles;
    for jj = 1:length(initAlignedParticlesRot)
        pNow = initAlignedParticlesRot{jj}.points;
        pNow = rmat*pNow';
        initAlignedParticlesRot{jj}.points = pNow';
    end

    drawnow

    %% Run the bootstrapping

%     try
% 
%         load([dataDirNow saveTag '_piezoAveragingWorkspaceRotated.mat'])
% 
%         if ~exist('initAlignedParticlesRot','var')
%             error('Loading failed')
%         end
% 
%     catch

        % bootstrapping with imposing symmetry prior knowledge
        USE_SYMMETRY = 0;   % flag for imposing symmetry prio knowledge (n-fold = #)
        M1 = [];            % not implemented
        [superParticleWithoutPKRot, ~] = one2all3D(initAlignedParticlesRot, iter, M1, [dataDir 'withSym2_' saveTag '_'], sup, USE_SYMMETRY, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);
    
        visualizeSMLM3D(superParticleWithoutPKRot{1,6},scale/2, 1);
        xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize',16)
        set(gcf, 'InvertHardcopy', 'off')
        saveas(gcf,[dataDir 'superParticleWithSym2_' saveTag '.png'],'png')
        saveas(gcf,[dataDir 'superParticleWithSym2_' saveTag '.fig'],'fig')
    
    
        % bootstrapping with imposing symmetry prior knowledge
        USE_SYMMETRY = 2;   % flag for imposing symmetry prio knowledge (n-fold = #)
        M1 = [];            % not implemented
        [superParticleWithPK2Rot, ~] = one2all3D(initAlignedParticlesRot, iter, M1, [dataDir 'withSym2_' saveTag '_'], sup, USE_SYMMETRY, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);
    
        visualizeSMLM3D(superParticleWithPK2Rot{1,6},scale/2, 1);
        xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize',16)
        set(gcf, 'InvertHardcopy', 'off')
        saveas(gcf,[dataDir 'superParticleWithSym2_' saveTag '.png'],'png')
        saveas(gcf,[dataDir 'superParticleWithSym2_' saveTag '.fig'],'fig')
    
        % bootstrapping with imposing symmetry prior knowledge
        USE_SYMMETRY = 3;   % flag for imposing symmetry prio knowledge (n-fold = #)
        M1 = [];            % not implemented
        [superParticleWithPK3Rot, ~] = one2all3D(initAlignedParticlesRot, iter, M1, [dataDir 'withSym4_' saveTag '_'], sup, USE_SYMMETRY, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);
    
        visualizeSMLM3D(superParticleWithPK3Rot{1,6},scale/2, 1);
        xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize',16)
        set(gcf, 'InvertHardcopy', 'off')
        saveas(gcf,[dataDir 'superParticleWithSym4_' saveTag '.png'],'png')
        saveas(gcf,[dataDir 'superParticleWithSym4_' saveTag '.fig'],'fig')
    
        % bootstrapping with imposing symmetry prior knowledge
        USE_SYMMETRY = 4;   % flag for imposing symmetry prio knowledge (n-fold = #)
        M1 = [];            % not implemented
        [superParticleWithPK4Rot, ~] = one2all3D(initAlignedParticlesRot, iter, M1, [dataDir 'withSym4_' saveTag '_'], sup, USE_SYMMETRY, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);
    
        visualizeSMLM3D(superParticleWithPK4Rot{1,6},scale/2, 1);
        xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize',16)
        set(gcf, 'InvertHardcopy', 'off')
        saveas(gcf,[dataDir 'superParticleWithSym4_' saveTag '.png'],'png')
        saveas(gcf,[dataDir 'superParticleWithSym4_' saveTag '.fig'],'fig')
    
        % bootstrapping with imposing symmetry prior knowledge
        USE_SYMMETRY = 5;   % flag for imposing symmetry prio knowledge (n-fold = #)
        M1 = [];            % not implemented
        [superParticleWithPK5Rot, ~] = one2all3D(initAlignedParticlesRot, iter, M1, [dataDir 'withSym4_' saveTag '_'], sup, USE_SYMMETRY, USE_GPU_GAUSSTRANSFORM, USE_GPU_EXPDIST);
    
        visualizeSMLM3D(superParticleWithPK5Rot{1,6},scale/2, 1);
        xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize',16)
        set(gcf, 'InvertHardcopy', 'off')
        saveas(gcf,[dataDir 'superParticleWithSym5_' saveTag '.png'],'png')
        saveas(gcf,[dataDir 'superParticleWithSym5_' saveTag '.fig'],'fig')

%     end

    %% Animate Results
    if mm <= 11
        % Manually set the axis to a range I like
        axLim = 40*[-1 1 -1 1 -1 1];
        cLim = [10^-5 8*10^-5];
    else
        axLim = [];
        cLim = [];
        saveTag = saveTag(10:end); % There is a file name length problem with these files
    end

    data = superParticleWithoutPKRot{1,6};
    figName = [saveTag '_noSymRot'];
    animateParticle(data,scale,figName,dataDirNow,axLim,cLim)

    data = superParticleWithoutPK{1,6};
    figName = [saveTag '_noSymNoRot'];
    animateParticle(data,scale,figName,dataDirNow,axLim,cLim)

    data = superParticleWithPK2Rot{1,6};
    figName = [saveTag '_2foldRot'];
    animateParticle(data,scale,figName,dataDirNow,axLim,cLim)

    data = superParticleWithPK3Rot{1,6};
    figName = [saveTag '_3foldRot'];
    animateParticle(data,scale,figName,dataDirNow,axLim,cLim)

    data = superParticleWithPK4Rot{1,6};
    figName = [saveTag '_4foldRot'];
    animateParticle(data,scale,figName,dataDirNow,axLim,cLim)

    data = superParticleWithPK5Rot{1,6};
    figName = [saveTag '_5foldRot'];
    animateParticle(data,scale,figName,dataDirNow,axLim,cLim)   

    %% Save everything
    save([dataDirNow saveTag '_piezoAveragingWorkspaceRotated.mat'])

end

%% WRAP UP
disp(' ')
disp('BATCH SCRIPT COMPLETE')