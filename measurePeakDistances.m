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
    
    h = histogram(densityGrid(:));
    T2 = otsuthresh(h.BinCounts);    
    T2 = T2*max(densityGrid(:));

    % Find pixels that are above the threshold
    k = find(densityGrid>T2);
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

    pkC = NaN*pkx;
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

    dPeak = mvksdensity(data,pkC,'Bandwidth',scale/2); %  This is the probability density ath the peaks
    dPeakCoarse = mvksdensity(data,pkC,'Bandwidth',scale/2); %  This is the probability density ath the peaks

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

    %% Save the data    
    save([directory saveTag '_PeakFinding.mat'],'data','density','blades','x','y','z','d',...
        'densityGrid','pkx','pkC','toPlotCoarse','toPlotFine','dPeak','dPeakCoarse')

end

%% Clean up
close all
disp('Batch Complete')