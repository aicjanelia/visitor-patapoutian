%%%%%% INSTRUCTIONS
% (1) Run batch_piezoSymFold over the data you wish to analyze. The
%       relevant file created by that process is
%       '*_piezoAveragingWorkspaceRotated.mat'
% (2) Fill in the user parameters below and run the script. The variable
%       directories allows this script to loop through multiple experiments,
%       analyzing each one in turn.

clc, clear, close all

%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directories = {'Y:\Rachel\Patapoutian_iPALM_data_analyzed\iPALM_data_processed_EM\CombinedParticles\negYODA1\',...
    'Y:\Rachel\Patapoutian_iPALM_data_analyzed\iPALM_data_processed_EM\CombinedParticles\posYODA1_2\'};

scale = 5; % This is the scale used in the averaging process
Nresample = 2000;
epsilon = 1; % This will jump over the pixel-lockign gaps
minpts = Nresample/4; % Higher density regions = shows up in more bootstrap instances, try 25%
CIalpha = 0.025; % 0.025*2 = 0.05 for 95% CI

%% MAKE MEASUREMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for mm = 1:length(directories)

    directory = directories{mm};
    fileName = dir([directory '*_piezoAveragingWorkspaceRotated.mat']);
    if isempty(fileName)
        error('Rotated SuperParticle Not Found')
    end
    fileName = fileName(1).name;
    saveTag = fileName(1:end-35);

    %% Load the raw data
    
    load([directory fileName])

    data = superParticleWithPK3Rot{1,6};    
    density = mvksdensity(data,data,'Bandwidth',scale/2); %  This is the probability density
    
    %% Calculate a threshold
    
    h = histogram(density(:));
    T = otsuthresh(h.BinCounts);
    blades = density > T*max(density(:));

    x = data(blades,1);
    y = data(blades,2);
    z = data(blades,3);
    d = density(blades);
    
    figure(1)
    set(gcf,'Position',[500 275 560*2 420*2])    
    scatter3(x,y,z, d*5*10^6, d, '.')
    colormap(cool)
    xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'FontSize',16)
    colorbar       
    set(gcf,'Color','white')
    
    view([46 49])
    xlim(25*[-1 1])
    ylim(25*[-1 1])
    zlim(10*[-1 1])
    caxis((10^-5)*[2 6])
    
    saveas(gcf,[directory saveTag '_OtsuSuperParticle.png'],'png')
    saveas(gcf,[directory saveTag '_OtsuSuperParticle.fig'],'fig')
    
    %% Create Gridded Data to Calculate Coarse Peaks
    
    meshPos = ceil(max([max(x) abs(min(x)) max(y) abs(min(y))]))+5; % Max position on the grid
    meshSpace = 1; % 1 is so indices = nm (could consider other options if necessary)
    
    [Xq,Yq,Zq] = meshgrid(-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos);
    Vq = mvksdensity(data,[Xq(:) Yq(:) Zq(:)],'Bandwidth',scale/2);    
    densityGrid = reshape(Vq,size(Xq));

    figure(2)
    imagesc(densityGrid(:,:,meshPos))
    set(gca,'DataAspectRatio',[1 1 1])
    saveas(gcf,[directory saveTag '_densityGrid_slice' num2str(meshPos) '.png'],'png')
    saveas(gcf,[directory saveTag '_densityGrid_slice' num2str(meshPos) '.fig'],'fig')

    % Find pixels that are above the threshold
    k = find(densityGrid>T*max(density(:)));
    [a,b,c]=ind2sub(size(densityGrid),k);
    edgeSize = max(size(densityGrid));

    keepers = 0*a; % Values that count as peaks
    for kk = 1:length(a)
        x1=a(kk);
        x2=b(kk);
        x3=c(kk);
        % This is an adaptation of the algorithm in pkfnd, adapted to 3D data
        if x1>1 && x1<edgeSize && x2>1 && x2<edgeSize && x3>1 && x3<edgeSize % avoid the image edges        
            if densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2,x3) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2,x3) && ... % all 26 neighbors
               densityGrid(x1,x2,x3) >= densityGrid(x1,x2-1,x3) && densityGrid(x1,x2,x3) >= densityGrid(x1,x2+1,x3) && ...    
               densityGrid(x1,x2,x3) >= densityGrid(x1,x2,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1,x2,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2-1,x3) && densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2+1,x3) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2-1,x3) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2+1,x3) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1,x2-1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1,x2-1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1,x2+1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1,x2+1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2-1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2-1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2+1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2+1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2-1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2-1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2+1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2+1,x3+1)
    
                    keepers(kk) = 1; % If the pixel is the brightest of all its neighbors, keep it
    
            end
    
        end
    end
    
    pkx = [a(keepers>0) b(keepers>0) c(keepers>0)]; % Peaks in the mesh grid coordinate system
    
    %% Calculate Fine Peaks Based on Point Cloud 

    sz = scale; %window for fitting neighbor peaks

    pkC = NaN*ones(size(pkx,1),3);
    for kk = 1:sum(keepers)

        r = sqrt(sum((data(:,1)-pkx(kk,2)+meshPos).^2 + (data(:,2)-pkx(kk,1)+meshPos).^2 + (data(:,3)-pkx(kk,3)+meshPos).^2,2));
        neigh = (r<sz);
        if ~sum(neigh)
            error('No neighbors found around peak')
        end

        tmp = density(neigh);
        norm = sum(tmp);
        xavg = sum(data(neigh,1).*tmp)./norm; % Weighted average based on local density
        yavg = sum(data(neigh,2).*tmp)./norm;
        zavg = sum(data(neigh,3).*tmp)./norm;

        pkC(kk,:) = [xavg yavg zavg]; % This is in nm relative in the superparticle coordinate system

    end
%     disp(pkC)

    %% Plot locations of found peaks
    
    toPlotCoarse = [b(keepers>0) a(keepers>0) c(keepers>0)]-meshPos; % Move back to particle coordinate system, row = y, col = x
    toPlotFine = [pkC(:,1) pkC(:,2) pkC(:,3)];
    
    figure(3)
    set(gcf,'Position',[500 275 560*2 420*2])    
    scatter3(x,y,z, d*5*10^6, d, '.')
    hold on
    l1 = plot3(toPlotCoarse(:,1),toPlotCoarse(:,2),toPlotCoarse(:,3),'ok','MarkerSize',10,'LineWidth',2);
    l2 = plot3(toPlotFine(:,1),toPlotFine(:,2),toPlotFine(:,3),'k^','MarkerSize',10,'LineWidth',2);
    hold off
    colormap(cool)
    xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'FontSize',16)
    colorbar       
    set(gcf,'Color','white')
    legend([l1 l2],'Coarse Peaks','Fine Peaks')
    
    view(2)
    
    saveas(gcf,[directory saveTag '_FoundPeaksOnParticle.png'],'png')
    saveas(gcf,[directory saveTag '_FoundPeaksOnParticle.fig'],'fig')

    dPeak = mvksdensity(data,pkC(:,1:3),'Bandwidth',scale/2); %  This is the probability density ath the peaks
    dPeakCoarse = mvksdensity(data,pkC(:,1:3),'Bandwidth',scale/2); %  This is the probability density ath the peaks

    figure(4)
    scatter3(data(:,1),data(:,2),density(:),[],density(:))
    set(gca,'DataAspectRatio',[1 1 10^-6])
    hold on
    for kk = 1:size(pkC,1)
        l1 = plot3(toPlotCoarse(kk,1),toPlotCoarse(kk,2),dPeakCoarse(kk),'ok','MarkerSize',10);
        l2 = plot3(pkC(kk,1),pkC(kk,2),dPeak(kk),'pr','MarkerFaceColor','r','MarkerSize',10);
    end
    hold off
    legend([l1 l2],'Coarse Peaks','Fine Peaks')
    saveas(gcf,[directory saveTag '_FoundPeaks2D.png'],'png')
    saveas(gcf,[directory saveTag '_FoundPeaks2D.fig'],'fig')

    %% Save the direct data as an intermediate
    save([directory saveTag '_PeakFinding.mat'],'data','density','blades','x','y','z','d',...
        'densityGrid','pkx','pkC','toPlotCoarse','toPlotFine','dPeak','dPeakCoarse')


    %% Intermediate clean up
    close all

    %% Bootstrap Sampling
    disp('Starting Bootstrap')

    s = RandStream('mlfg6331_64'); % Reproducible random stream
    indSample = 1:size(data,1);
    pkCBoot = cell(Nresample,1);
    dataBoot = pkCBoot; densityBoot = pkCBoot; bladesBoot = pkCBoot;
    xBoot = pkCBoot; yBoot = xBoot; zBoot = xBoot; dBoot = xBoot;
    dPeakBoot = pkCBoot;

    for ii = 1:Nresample
        if ~mod(ii,50)
            disp(ii)
        end

        %%% Randomly sample with replacement to a dataset of the same size
        dataSamp = datasample(s,indSample,size(data,1),'Replace',true);
        dataSamp = data(dataSamp,:);
        
        densitySamp = mvksdensity(dataSamp,dataSamp,'Bandwidth',scale/2); %  This is the probability density
    
        %%% Calculate a threshold        
        figure(10)
        h = histogram(densitySamp(:));
        T = otsuthresh(h.BinCounts);
        bladesSamp = densitySamp > T*max(densitySamp(:));
    
        xSamp = dataSamp(bladesSamp,1);
        ySamp = dataSamp(bladesSamp,2);
        zSamp = dataSamp(bladesSamp,3);
        dSamp = densitySamp(bladesSamp);
        
%         figure(1)
%         set(gcf,'Position',[500 275 560*2 420*2])    
%         scatter3(xSamp,ySamp,zSamp, dSamp*5*10^6, dSamp, '.')
%         colormap(cool)
%         xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
%         set(gca,'DataAspectRatio',[1 1 1])
%         set(gca,'FontSize',16)
%         colorbar       
%         set(gcf,'Color','white')
%         drawnow
%         
%         view([46 49])
%         xlim(25*[-1 1])
%         ylim(25*[-1 1])
%         zlim(10*[-1 1])
%         caxis((10^-5)*[2 6])
    
        %%% Create Gridded Data to Calculate Coarse Peaks        
        meshPos = ceil(max([max(xSamp) abs(min(xSamp)) max(ySamp) abs(min(ySamp))]))+5; % Max position on the grid
        meshSpace = 1; % 1 is so indices = nm (could consider other options if necessary)
        
        [Xq,Yq,Zq] = meshgrid(-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos);
        Vq = mvksdensity(dataSamp,[Xq(:) Yq(:) Zq(:)],'Bandwidth',scale/2);    
        densityGridSamp = reshape(Vq,size(Xq));
    
%         figure(2)
%         imagesc(densityGridSamp(:,:,meshPos))
%         set(gca,'DataAspectRatio',[1 1 1])

        % Don't allow isolated components, only want to keep the main PIEZO
        % structure, not small peaks out in the 'noise'
        vol = densityGridSamp>T*max(density(:));
        CC = bwconncomp(vol);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        filtered_vol = false(size(vol));
        filtered_vol(CC.PixelIdxList{idx}) = true;
    
        % Find pixels that are above the threshold
        k = find(filtered_vol);
        [a,b,c]=ind2sub(size(densityGridSamp),k);
        edgeSize = max(size(densityGridSamp));
    
        keepers = 0*a; % Values that count as peaks
        for kk = 1:length(a)
            x1=a(kk);
            x2=b(kk);
            x3=c(kk);
            % This is an adaptation of the algorithm in pkfnd, adapted to 3D data
            if x1>1 && x1<edgeSize && x2>1 && x2<edgeSize && x3>1 && x3<edgeSize % avoid the image edges        
                if densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3) && ... % all 26 neighbors
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3) && ...    
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3+1)
        
                        keepers(kk) = 1; % If the pixel is the brightest of all its neighbors, keep it
        
                end
        
            end
        end
        
        pkxSamp = [a(keepers>0) b(keepers>0) c(keepers>0)]; % Peaks in the mesh grid coordinate system
        
        %%% Calculate Fine Peaks Based on Point Cloud     
        sz = scale; %window for fitting neighbor peaks
    
        pkCSamp = NaN*ones(size(pkxSamp,1),3);
        for kk = 1:sum(keepers)
    
            r = sqrt(sum((dataSamp(:,1)-pkxSamp(kk,2)+meshPos).^2 + (dataSamp(:,2)-pkxSamp(kk,1)+meshPos).^2 + (dataSamp(:,3)-pkxSamp(kk,3)+meshPos).^2,2));
            neigh = (r<sz);
            if ~sum(neigh)
                error('No neighbors found around peak')
            end
    
            tmp = densitySamp(neigh);
            norm = sum(tmp);
            xavg = sum(dataSamp(neigh,1).*tmp)./norm; % Weighted average based on local density
            yavg = sum(dataSamp(neigh,2).*tmp)./norm;
            zavg = sum(dataSamp(neigh,3).*tmp)./norm;
    
            pkCSamp(kk,:) = [xavg yavg zavg]; % This is in nm relative in the superparticle coordinate system
    
        end

        figure(3)
        set(gcf,'Position',[500 275 560*2 420*2])    
        scatter3(xSamp,ySamp,zSamp, dSamp*5*10^6, dSamp, '.')
        hold on
        plot3(pkCSamp(:,1),pkCSamp(:,2),pkCSamp(:,3),'pr','MarkerFaceColor','r','MarkerSize',10);
        hold off
        colormap(cool)
        xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize',16)
        colorbar       
        set(gcf,'Color','white')
        view(2)
        drawnow

        dPeakSamp = mvksdensity(dataSamp,pkCSamp(:,1:3),'Bandwidth',scale/2); %  This is the probability density at the peaks
        figure(4)
        scatter3(dataSamp(:,1),dataSamp(:,2),densitySamp(:),[],densitySamp(:))
        set(gca,'DataAspectRatio',[1 1 10^-6])
        hold on
        for kk = 1:size(pkCSamp,1)
            plot3(pkCSamp(kk,1),pkCSamp(kk,2),dPeakSamp(kk),'pr','MarkerFaceColor','r','MarkerSize',10);
        end
        hold off
        drawnow

        pkCBoot{ii} = pkCSamp;

        dataBoot{ii} = dataSamp;
        densityBoot{ii} = densitySamp;
        bladesBoot{ii} = bladesSamp;
        xBoot{ii} = xSamp;
        yBoot{ii} = ySamp;
        zBoot{ii} = zSamp;
        dBoot{ii} = dSamp;
        dPeakBoot{ii} = dPeakSamp;

    end

    save([directory saveTag '_PeakFindingBootstrap.mat'],'data','density','blades','x','y','z','d',...
        'densityGrid','pkx','pkC','toPlotCoarse','toPlotFine','dPeak','dPeakCoarse',...
        'pkCBoot','dataBoot','densityBoot','bladesBoot','xBoot','yBoot','zBoot','dBoot','dPeakBoot')

    %% Leave One Out Analysis...    
    disp('Starting LOO')

    Nloo = size(data,1);
    pkCLoo = cell(Nloo,1);
    dataLoo = pkCLoo; densityLoo = pkCLoo; bladesLoo = pkCLoo;
    xLoo = pkCLoo; yLoo = xLoo; zLoo = xLoo; dLoo = xLoo;
    dPeakLoo = pkCLoo;

    for ii = 1:Nloo
        if ~mod(ii,50)
            disp(ii)
        end

        % Leave one out
        dataSamp = ones(Nloo,1);
        dataSamp(ii) = 0;
        dataSamp = data(dataSamp>0,:);
        
        densitySamp = mvksdensity(dataSamp,dataSamp,'Bandwidth',scale/2); %  This is the probability density
    
        %%% Calculate a threshold        
        figure(10)
        h = histogram(densitySamp(:));
        T = otsuthresh(h.BinCounts);
        bladesSamp = densitySamp > T*max(densitySamp(:));
    
        xSamp = dataSamp(bladesSamp,1);
        ySamp = dataSamp(bladesSamp,2);
        zSamp = dataSamp(bladesSamp,3);
        dSamp = densitySamp(bladesSamp);
        
%         figure(1)
%         set(gcf,'Position',[500 275 560*2 420*2])    
%         scatter3(xSamp,ySamp,zSamp, dSamp*5*10^6, dSamp, '.')
%         colormap(cool)
%         xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
%         set(gca,'DataAspectRatio',[1 1 1])
%         set(gca,'FontSize',16)
%         colorbar       
%         set(gcf,'Color','white')
%         drawnow
%         
%         view([46 49])
%         xlim(25*[-1 1])
%         ylim(25*[-1 1])
%         zlim(10*[-1 1])
%         caxis((10^-5)*[2 6])
    
        %%% Create Gridded Data to Calculate Coarse Peaks        
        meshPos = ceil(max([max(xSamp) abs(min(xSamp)) max(ySamp) abs(min(ySamp))]))+5; % Max position on the grid
        meshSpace = 1; % 1 is so indices = nm (could consider other options if necessary)
        
        [Xq,Yq,Zq] = meshgrid(-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos);
        Vq = mvksdensity(dataSamp,[Xq(:) Yq(:) Zq(:)],'Bandwidth',scale/2);    
        densityGridSamp = reshape(Vq,size(Xq));
    
%         figure(2)
%         imagesc(densityGridSamp(:,:,meshPos))
%         set(gca,'DataAspectRatio',[1 1 1])

        % Don't allow isolated components, only want to keep the main PIEZO
        % structure, not small peaks out in the 'noise'
        vol = densityGridSamp>T*max(density(:));
        CC = bwconncomp(vol);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        filtered_vol = false(size(vol));
        filtered_vol(CC.PixelIdxList{idx}) = true;
    
        % Find pixels that are above the threshold
        k = find(filtered_vol);
        [a,b,c]=ind2sub(size(densityGridSamp),k);
        edgeSize = max(size(densityGridSamp));
    
        keepers = 0*a; % Values that count as peaks
        for kk = 1:length(a)
            x1=a(kk);
            x2=b(kk);
            x3=c(kk);
            % This is an adaptation of the algorithm in pkfnd, adapted to 3D data
            if x1>1 && x1<edgeSize && x2>1 && x2<edgeSize && x3>1 && x3<edgeSize % avoid the image edges        
                if densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3) && ... % all 26 neighbors
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3) && ...    
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3+1)
        
                        keepers(kk) = 1; % If the pixel is the brightest of all its neighbors, keep it
        
                end
        
            end
        end
        
        pkxSamp = [a(keepers>0) b(keepers>0) c(keepers>0)]; % Peaks in the mesh grid coordinate system
        
        %%% Calculate Fine Peaks Based on Point Cloud 
    
        sz = scale; %window for fitting neighbor peaks
    
        pkCSamp = NaN*ones(size(pkxSamp,1),3);
        for kk = 1:sum(keepers)
    
            r = sqrt(sum((dataSamp(:,1)-pkxSamp(kk,2)+meshPos).^2 + (dataSamp(:,2)-pkxSamp(kk,1)+meshPos).^2 + (dataSamp(:,3)-pkxSamp(kk,3)+meshPos).^2,2));
            neigh = (r<sz);
            if ~sum(neigh)
                error('No neighbors found around peak')
            end
    
            tmp = densitySamp(neigh);
            norm = sum(tmp);
            xavg = sum(dataSamp(neigh,1).*tmp)./norm; % Weighted average based on local density
            yavg = sum(dataSamp(neigh,2).*tmp)./norm;
            zavg = sum(dataSamp(neigh,3).*tmp)./norm;
    
            pkCSamp(kk,:) = [xavg yavg zavg]; % This is in nm relative in the superparticle coordinate system
    
        end

%         figure(3)
%         set(gcf,'Position',[500 275 560*2 420*2])    
%         scatter3(xSamp,ySamp,zSamp, dSamp*5*10^6, dSamp, '.')
%         hold on
%         plot3(pkCSamp(:,1),pkCSamp(:,2),pkCSamp(:,3),'pr','MarkerFaceColor','r','MarkerSize',10);
%         hold off
%         colormap(cool)
%         xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
%         set(gca,'DataAspectRatio',[1 1 1])
%         set(gca,'FontSize',16)
%         colorbar       
%         set(gcf,'Color','white')
%         view(2)
%         drawnow
% 
%         dPeakSamp = mvksdensity(dataSamp,pkCSamp(:,1:3),'Bandwidth',scale/2); %  This is the probability density at the peaks
%         figure(4)
%         scatter3(dataSamp(:,1),dataSamp(:,2),densitySamp(:),[],densitySamp(:))
%         set(gca,'DataAspectRatio',[1 1 10^-6])
%         hold on
%         for kk = 1:size(pkCSamp,1)
%             plot3(pkCSamp(kk,1),pkCSamp(kk,2),dPeakSamp(kk),'pr','MarkerFaceColor','r','MarkerSize',10);
%         end
%         hold off
%         drawnow

        pkCLoo{ii} = pkCSamp;

        dataLoo{ii} = dataSamp;
        densityLoo{ii} = densitySamp;
        bladesLoo{ii} = bladesSamp;
        xLoo{ii} = xSamp;
        yLoo{ii} = ySamp;
        zLoo{ii} = zSamp;
        dLoo{ii} = dSamp;
        dPeakLoo{ii} = dPeakSamp;

    end

    save([directory saveTag '_PeakFindingLOO.mat'],'data','density','blades','x','y','z','d',...
        'densityGrid','pkx','pkC','toPlotCoarse','toPlotFine','dPeak','dPeakCoarse',...
        'pkCLoo','dataLoo','densityLoo','bladesLoo','xLoo','yLoo','zLoo','dLoo','dPeakLoo')


    %% Cluster and Connect to Original Peaks

    % Leave one out
    looAll = cell2mat(pkCLoo);
    looID = kmeans(looAll,size(pkC,1)); % Jackknifing should be very similar to the original based on high n, keep the same number of clusters as there were peaks
    looCent = NaN*ones(max(looID),3); % centroid of the cluster
    for cc = 1:max(looID(:))
        rowsNow = (looID == cc);
        looCent(cc,:) = mean(looAll(rowsNow,1:3),'omitnan');
    end

    % Bootstrap
    allPeaks = cell2mat(pkCBoot);
    allH = cell2mat(dPeakBoot);
    idClus = dbscan(allPeaks(:,1:3),epsilon,minpts);
    if max(idClus(:)) > size(pkC,1)
        error('Unexpected Mapping')
    end
    bootCent = NaN*ones(max(idClus),3);
    for cc = 1:max(idClus)
        rowsNow = idClus == cc;
        bootCent(cc,:) = mean(allPeaks(rowsNow,1:3));
    end

    % Connect clusters
    orig2clus = NaN*ones(size(pkC,1),1);
    orig2loo = orig2clus;
    for cc = 1:size(pkC,1)
        rowsNow = pkC(cc,1:3);
        orig2clus(cc) = knnsearch(bootCent,rowsNow);
        orig2loo(cc) = knnsearch(looCent,rowsNow);
    end
    clus2orig = NaN*ones(size(pkC,1),1);
    for cc = 1:size(pkC,1)
        rowsNow = bootCent(cc,1:3);
        clus2orig(cc) = knnsearch(pkC,rowsNow);
    end
    loo2orig = NaN*ones(size(pkC,1),1);
    for cc = 1:size(pkC,1)
        rowsNow = looCent(cc,1:3);
        loo2orig(cc) = knnsearch(pkC,rowsNow);
    end

    %% Plot Clustering Results

    figure(10)
    set(gcf,'Position',[400 400 560*2 420*2])
    scatter3(allPeaks(:,1),allPeaks(:,2),allPeaks(:,3),6,allH,'filled')
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)
    h = colorbar;
    set(get(h,'Label'),'String','Peak Height')    
    saveas(gcf,[directory saveTag 'BootstrapAllPeaks.png'],'png')

    figure(11)    
    set(gcf,'Position',[800 300 560*2 420*2])
    scatter3(allPeaks(:,1),allPeaks(:,2),allPeaks(:,3),6,idClus,'filled')
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)
    colormap([0.75 0.75 0.75 ; 0.75 0.75 0.75 ; lines(max(idClus(:)))])
    caxis([-1 max(idClus(:))])
    h = colorbar;
    set(get(h,'Label'),'String','DBSCAN cluster')    
    saveas(gcf,[directory saveTag 'BootstrapDBSCAN.png'],'png')

    figure(12)
    set(gcf,'Position',[1200 200 560*2 420*2])
    scatter3(looAll(:,1),looAll(:,2),looAll(:,3),[],looID)
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)
    colormap(lines(max(looID(:))))
    caxis([1 max(looID(:))])
    h = colorbar;
    set(get(h,'Label'),'String','kmeans cluster')    
    saveas(gcf,[directory saveTag 'LOOkmeans.png'],'png')

    %%   Uncertainty in Peak Locations
    alphaPeaks = NaN*ones(2,3,max(idClus));
    CIPeaks = alphaPeaks;
    ahatPeaks = NaN*ones(max(idClus),3);
    z0Peaks = ahatPeaks;
    for cc = 1:size(pkC,1)

        % First calculate the acceleration term
        sliceLoo = (looID == orig2loo(cc));
        sliceLoo = looAll(sliceLoo,1:3);
        Ui = (size(sliceLoo,1)-1)*(mean(sliceLoo)-sliceLoo);
        ahatPeaks(cc,:) = sum(Ui.^3)./(sum(Ui.^2).^3/2)/6;

        % Now calculate the bias term    
        rowsNow = (idClus == orig2clus(cc));
        z0Peaks(cc,:) = norminv(sum(allPeaks(rowsNow,1:3)<pkC(cc,1:3))/sum(rowsNow));

        % Combine to find the BCa alphas
        clear alpha
        alpha(1,:) = normcdf( z0Peaks(cc,:) + (z0Peaks(cc,:) + norminv(CIalpha)) / (1 - ahatPeaks(cc,:).*(z0Peaks(cc,:) + norminv(CIalpha))) );
        alpha(2,:) = normcdf( z0Peaks(cc,:) + (z0Peaks(cc,:) + norminv(1-CIalpha)) / (1 - ahatPeaks(cc,:).*(z0Peaks(cc,:) + norminv(1-CIalpha))) );
        alphaPeaks(:,:,cc) = alpha;

        % Use BCa alphas to construct CI    
        CI = NaN*alpha;
        for dd = 1:size(alpha,2)
            CI(:,dd) = prctile(allPeaks(rowsNow,dd),alpha(:,dd)*100);
        end
        CIPeaks(:,:,cc) = CI;

        figure(cc+100)
        set(gcf,'Position',[200 200 560*2.5 420])
        subplot(1,3,1)
        histogram(allPeaks(rowsNow,1));
        hold on
        plot(pkC(cc,1)*[1 1],[0 300],'Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(1,1)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(2,1)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(mean(allPeaks(rowsNow,1))*[1 1],[0 300],'LineWidth',2)
        hold off        
        title('x')

        subplot(1,3,2)
        histogram(allPeaks(rowsNow,2));
        hold on
        plot(pkC(cc,2)*[1 1],[0 300],'Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(1,2)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(2,2)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(mean(allPeaks(rowsNow,2))*[1 1],[0 300],'LineWidth',2)
        hold off
        title('y')

        subplot(1,3,3)
        histogram(allPeaks(rowsNow,3));
        hold on
        plot(pkC(cc,3)*[1 1],[0 300],'Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(1,3)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(2,3)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(mean(allPeaks(rowsNow,3))*[1 1],[0 300],'LineWidth',2)
        hold off
        title('z')

        saveas(gcf,[directory saveTag 'BootstrapDist_Peak' num2str(cc) '.png'],'png')

    end

    % Figure of the peaks with confidence intervals
    figure(13)
    set(gcf,'Position',[1200 200 560*2 420*2])
    plot3(pkC(:,1),pkC(:,2),pkC(:,3),'ok','MarkerFaceColor','k','MarkerSize',6)
    hold on
    for cc = 1:size(CIPeaks,3)
        plot3(CIPeaks(:,1,cc),pkC(cc,2)*[1 1],pkC(cc,3)*[1 1],'k','LineWidth',2)
        plot3(pkC(cc,1)*[1 1],CIPeaks(:,2,cc),pkC(cc,3)*[1 1],'k','LineWidth',2)
        plot3(pkC(cc,1)*[1 1],pkC(cc,2)*[1 1],CIPeaks(:,3,cc),'k','LineWidth',2)
    end    
    hold off
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)
    saveas(gcf,[directory saveTag 'peakLocationCI.png'],'png')


    %% Uncertainty in Distances

    % For each instance, want to calculate the distances, so need to keep
    % track of which instance each peak come from
    eachPeak = NaN*ones(size(allPeaks,1),5);
    counter = 1;    
    for cc = 1:length(pkCBoot)
        slice = pkCBoot{cc}(:,1:3);
        eachPeak(counter:counter+size(slice,1)-1,:) = [slice cc*ones(size(slice,1),1) dPeakBoot{cc}];
        counter = counter+size(slice,1);
    end

    nPeaks = size(pkC,1);
    % Boostrap Distances  
    distBoot = NaN*ones(nPeaks,nPeaks,length(pkCBoot));
    for cc = 1:length(pkCBoot)
    
        slice = eachPeak(:,end-1) == cc;
        idSlice = idClus(slice);
        slice = eachPeak(slice,:);
    
        % Remove the noise cluster from consideration
        noiseClus = (idSlice == -1);
        slice = slice(~noiseClus,:);
        idSlice = idSlice(~noiseClus);
        idSlice = clus2orig(idSlice); % move into the pkC reference
    
        toPdist = NaN*ones(nPeaks,3);
        for dd = 1:nPeaks
            sliceNow = idSlice == dd;
            if sum(sliceNow) > 1
                sliceNow = slice(sliceNow,1:3);
                warning(['Extra Peaks detected for instance ' num2str(cc) ' (' num2str(size(sliceNow,1)) ' on peak ' num2str(dd) ')'])
                toPdist(dd,:) = mean(sliceNow); % TODO: Is this the right behavior for a duplicate?
            elseif sum(sliceNow)
                toPdist(dd,:) = slice(sliceNow,1:3);
            end
        end
    
        distBoot(:,:,cc) = squareform(pdist(toPdist));
    
    end

    % For each LOO, want to calculate the distances, so need to keep
    % track of which instance each peak come from
    eachLoo = NaN*ones(size(looAll,1),5);
    counter = 1;    
    for cc = 1:length(pkCLoo)
        slice = pkCLoo{cc}(:,1:3);
        eachLoo(counter:counter+size(slice,1)-1,:) = [slice cc*ones(size(slice,1),1) dPeakLoo{cc}];
        counter = counter+size(slice,1);
    end

    % LOO Distances  
    distLoo = NaN*ones(nPeaks,nPeaks,length(pkCLoo));
    for cc = 1:length(pkCLoo)
    
        slice = eachLoo(:,end-1) == cc;
        idSlice = looID(slice);
        slice = eachLoo(slice,:);
    
        % Remove the noise cluster from consideration
        noiseClus = (idSlice == -1);
        slice = slice(~noiseClus,:);
        idSlice = idSlice(~noiseClus);
        idSlice = loo2orig(idSlice);
    
        toPdist = NaN*ones(nPeaks,3);
        for dd = 1:nPeaks
            sliceNow = idSlice == dd;
            if sum(sliceNow) > 1
                sliceNow = slice(sliceNow,1:3);
                warning(['Extra LOO Peaks detected for instance ' num2str(cc) ' (' num2str(sliceNow) ' on peak ' num2str(dd) ')'])
                toPdist(dd,:) = mean(sliceNow); % TODO: Is this the right behavior for a duplicate?
            elseif sum(sliceNow)
                toPdist(dd,:) = slice(sliceNow,1:3);
            end
        end
    
        distLoo(:,:,cc) = squareform(pdist(toPdist));
    
    end

    %%% Confidence Intervals
    alphaDist = NaN*ones(2,nPeaks,nPeaks);
    CIDist = alphaDist;
    ahatDist = NaN*ones(nPeaks);
    z0Dist = ahatDist;

    figure(21) % Will be all the relevant distances
    set(gcf,'Position',[400 400 560*2 420*2])
    for cc = 1:nPeaks
        for dd = cc+1:nPeaks

            dOrig = sqrt(sum((pkC(cc,1:3)-pkC(dd,1:3)).^2));

            % First calculate the acceleration term
            sliceLoo = squeeze(distLoo(cc,dd,:));
            sliceLoo = sliceLoo(~isnan(sliceLoo)); % In case of any instances where the peak was missing
            Ui = (size(sliceLoo,1)-1)*(mean(sliceLoo)-sliceLoo);
            ahatDist(cc,dd) = sum(Ui.^3)./(sum(Ui.^2).^3/2)/6;

            % Now calculate the bias term    
            sliceBoot = squeeze(distBoot(cc,dd,:));
            z0Dist(cc,dd) = norminv(sum(sliceBoot<dOrig)/sum(~isnan(sliceBoot)));

            % Combine to find the BCa alphas
            clear alpha
            alpha(1) = normcdf( z0Dist(cc,dd) + (z0Dist(cc,dd) + norminv(CIalpha)) / (1 - ahatDist(cc,dd).*(z0Dist(cc,dd) + norminv(CIalpha))) );
            alpha(2) = normcdf( z0Dist(cc,dd) + (z0Dist(cc,dd) + norminv(1-CIalpha)) / (1 - ahatDist(cc,dd).*(z0Dist(cc,dd) + norminv(1-CIalpha))) );
            alphaDist(:,cc,dd) = alpha;

            % Use BCa alphas to construct CI    
            CI = NaN*alpha;
            for ee = 1:length(alpha)
                CI(ee) = prctile(sliceBoot,alpha(:,ee)*100);
            end
            CIDist(:,cc,dd) = CI;

            subplot(nPeaks-1,nPeaks-1,(nPeaks-1)*(cc-1)+dd-1)
            histogram(distBoot(cc,dd,:),15)
            hold on
            plot(dOrig*[1 1],[0 300],'Color',[0.8500    0.3250    0.0980],'LineWidth',2)
            plot(CI(1)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
            plot(CI(2)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
            plot(mean(sliceBoot,'omitnan')*[1 1],[0 300],'LineWidth',2)
            hold off        
            title([num2str(cc) ' - ' num2str(dd)])
        end
    end
    saveas(gcf,[directory saveTag 'BootstrapDist_EachDistance.png'],'png')

    %% Figure of Distances
    % Figure of the peaks with confidence intervals
    dOrig = squareform(pdist(pkC));
    cmapPlot = lines(nPeaks.^2);

    figure(22)
    set(gcf,'Position',[400 200 560*2 420*2])

    plot3(pkC(:,1),pkC(:,2),pkC(:,3),'ok','MarkerFaceColor','k','MarkerSize',8)

    hold on
    for cc = 1:nPeaks
        for dd = cc+1:nPeaks
            plot3([pkC(cc,1) pkC(dd,1)],[pkC(cc,2) pkC(dd,2)],[pkC(cc,3) pkC(dd,3)],'Color',cmapPlot(nPeaks*(cc-1)+dd,:),'LineWidth',2)
            toPlot = mean(pkC([cc dd],:));
            text(toPlot(:,1),toPlot(:,2),toPlot(:,3),...
                [num2str(dOrig(cc,dd),'%0.1f') ' (' num2str(CIDist(1,cc,dd),'%0.1f') ',' num2str(CIDist(2,cc,dd),'%0.1f') ')'],...
                'Color',cmapPlot(nPeaks*(cc-1)+dd,:),'FontSize',16)
        end
    end    
    hold off

    xlabel('x','FontSize',16)
    ylabel('y','FontSize',16)
    zlabel('z','FontSize',16)
    set(gca,'FontSize',16)
    saveas(gcf,[directory saveTag 'distanceCI.fig'],'fig') % This figure makes no sense static, save fig for exploration   


    %% Save the data    
    save([directory saveTag '_PeakFindingAll.mat']) % This is likely to be an annoying large file
    
end


%% Clean up
% close all
disp('Batch Complete')