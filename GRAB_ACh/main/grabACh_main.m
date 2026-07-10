%% ========================================================================
%  1. User settings and paths
%  ========================================================================

% Get the path to this script and the repository root.
% Assumes this file is inside:
%   GRAB-ACh/main/
mainDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(mainDir);

% Add repository functions to the MATLAB path.
addpath(genpath(fullfile(repoRoot, 'functions')));

% Optional: add external dependencies if you use this folder.
externalDir = fullfile(repoRoot, 'external');
if exist(externalDir, 'dir')
    addpath(genpath(externalDir));
end

% -------------------------------------------------------------------------
% Session settings
% -------------------------------------------------------------------------

pat = 'Z:\Mellor Lab data\Dan Goodwin\Cohort 5\240307\DG21';
behaviourRoot = 'C:\Users\dgood\OneDrive - University of Exeter\Desktop\VR Data\Cohort 6\Behaviour';

animalID = 'DG21';
exptDate = '240307';
exptNo = 2;
fn = 'file_00001';
vrWithimaging = 1;
imgNum = 2;

% -------------------------------------------------------------------------
% Imaging settings
% -------------------------------------------------------------------------

Fs = 30;
dsTfactor = 1;
dsSfactor = 2;
targetDur = 5;

%% ========================================================================
%  2. Load raw ScanImage TIFF stack
%  ========================================================================
if ~strcmp(pat(end),'\'), pat = [pat '\']; end
pat = [pat num2str(exptNo) '\'];

disp('Reading tiff stack')   
%I = loadtiff([pat fn '.tif']);                     % alternative (slightly slower) method
% use multi-image tif using ScanImage tiff reader ---- % 
% see: https://vidriotech.gitlab.io/scanimagetiffreader-matlab/ 

% % % I = tiffreadVolume('D:\230915\DG12\1\file_00001.tif', 'PixelRegion', {[1 Inf], [1 Inf], [1 65000]}); % to import specific frames from the
% multi-image tiff

import ScanImageTiffReader.ScanImageTiffReader
reader = ScanImageTiffReader([pat fn '.tif']);      % create image object
I = reader.data();                                  % read frames from image object
[rPix,cPix,nFrames] = size(I);                      % size of image array
% A problem is that frames imported using the tiff reader are in the 
% wrong orientation and need to be transposed
for k = 1:nFrames                                   % cycle through each frame
    I(:,:,k) = flipud(rot90(I(:,:,k)));             % transpose to correct orientation and replace in original array
end

% % % I = uint16(I);

%% ========================================================================
%  3. Downsample image stack
%  ========================================================================
Fsds = Fs/dsTfactor;                                % downsampled sampling frequency (in Hz)
bin = 1:dsTfactor:nFrames-dsTfactor;                % bins for temporal downsampling
Ids = zeros(rPix/dsSfactor,cPix/dsSfactor,length(bin),'int16'); % pre-allocate downsampled output array
w0 = waitbar(0, 'Downsampling tiff stack');
for k = 1:length(bin)
    frames = I(:,:,bin(k):bin(k)+dsTfactor-1);      % frames in kth bin 
    framesAvg = mean(frames,3);                     % average time projection for bin frames
    framesAvgRs = imresize(framesAvg,1/dsSfactor);  % spatial downsampling via nearest neighbour interpolation 
    Ids(:,:,k) = int16(framesAvgRs);                % add downsampled image to output array
    waitbar(k/length(bin))
end 
close(w0)

% ---- Motion correction (subpixel DFT method (Thevenaz et al.) ---- %
if ~exist('targetDur','var'), targetDur = 1; end
target = mean(Ids(:,:,1:Fsds*targetDur), 3);        % target image for registration 

[rowPix, colPix, numFrames] = size(Ids);            % size of downsampled image array
precision = 100;                                    % subpixel precision for registration

w0 = waitbar(0, 'Registering tiff stack');
for k = 1:numFrames                                 % cyle through each frame
    frame = Ids(:,:,k);                             % kth frame 
    [transformOutput,tmpIreg] = dftregistration(fft2(target),fft2(frame),precision);  % register kth frame to target using DFT method
    Ireg = abs(ifft2(tmpIreg));                     % reconstruct registered frame k
    Ids(:,:,k) = Ireg;                              % add registered frame k to array containing all registered frames
    waitbar(k/numFrames)
end
close(w0)

trace = zeros(numFrames,1);                            % pre-allocate
for k = 1:numFrames                                    % cycle through frames
    frame = Ids(:,:,k);                        % kth frame
    trace(k) = mean(frame(:));                      % mean of kth frame
end

%% --- Create mask for GRAB-ACh signal (i.e. omitting blood vessels, etc) --- % - manually picking GRAB-ACh mask and background mask
Idsm = mean(Ids,3);                                 % mean projection 
prct = 20;                                          % upper N percent of signal used for GRAB-ACh region 
% thresh = min(Idsm(:)) + (max(Idsm(:))-min(Idsm(:)))*(prct/100); % image threshold 
[thresh,regions] = slider_fig_GRAB(Idsm);           % set threshold manually

Idsmbw = Idsm>thresh;                               % binarised image 

% ---- Approximate background signal (i.e. blood vessels) ---- % 
prctBg = 2;                                         % lowest N percent of pixels                               
threshBg = prctile(sort(Idsm(:),'ascend'),prctBg);
% IdsmBg = Idsm<threshBg;                             % binarised background image 
[IdsmBg,roiX,roiY] = drawROI_GRAB(Idsm);             % binarised background image 
%%
Fbg = zeros(numFrames,1);                           % pre-allocate raw fluorescence trace array
for k = 1:numFrames         
    frame = Ids(:,:,k);                             % kth image frame
    Fbg(k) = mean(frame(IdsmBg));                   % mean signal within background mask of frame k
end

% Fbg = zeros(numFrames,1); % quick way to remove the background subtraction from analysis - not convinced this is the best way to go yet

% calculate time vector for fluorescence trace  ---- %
dt = 1/Fsds;                                        % delta time between samples (in seconds)
tvec = dt:dt:dt*numFrames;                          % time vector

%%
% ---- Calculate deltaF/F GRAB-ACh fluorescence trace  ---- %
F = zeros(numFrames,1);                             % pre-allocate raw fluorescence trace array
for k = 1:numFrames         
    frame = Ids(:,:,k);                             % kth image frame
    F(k) = mean(frame(Idsmbw));                     % mean signal with GRAB-ACh mask area of frame k
end

FbgSub = F-Fbg; 

%% ========================================================================
%  6. Compare background-filtering approaches
%  ======================================================================== 

% Assume Fbg, F, numFrames, and tvec are already defined

% 1. Low-Pass Filter
cutoff_frequency = 1; % Adjust based on the expected frequencies in the background signal
Fbg_lowpass = lowpass(Fbg, cutoff_frequency, Fsds);

% 2. Savitzky-Golay Filter
polynomial_order = 1;  % Adjust based on the signal's smoothness
frame_size = 61;       % Window size for smoothing
Fbg_sgolay = sgolayfilt(Fbg, polynomial_order, frame_size);

% 3. Moving Average Filter (window size of 15 frames)
window_size = 61;
Fbg_moving_avg = movmean(Fbg, window_size);

% Plotting the background traces (original and filtered)
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]); % Wide figure

% Subplot 1: Original Background (Fbg)
subplot(4,1,1);
plot(tvec, Fbg, 'm', 'LineWidth', 1.5);
title('Original Background Signal (Fbg)');
ylabel('Fluorescence (A.U.)');
grid on;

% Subplot 2: Low-Pass Filtered Background
subplot(4,1,2);
plot(tvec, Fbg_lowpass, 'r', 'LineWidth', 1.5);
title('Low-Pass Filtered Background');
ylabel('Fluorescence (A.U.)');
grid on;

% Subplot 3: Savitzky-Golay Filtered Background
subplot(4,1,3);
plot(tvec, Fbg_sgolay, 'g', 'LineWidth', 1.5);
title('Savitzky-Golay Filtered Background');
ylabel('Fluorescence (A.U.)');
grid on;

% Subplot 4: Moving Average Filtered Background
subplot(4,1,4);
plot(tvec, Fbg_moving_avg, 'b', 'LineWidth', 1.5);
title('Moving Average Filtered Background');
xlabel('Time (s)');
ylabel('Fluorescence (A.U.)');
grid on;

% Optimize the layout
sgtitle('Comparison of Original and Filtered Background Signals');

% Save the background filter comparison figure
saveas(gcf, fullfile(pat, 'Background_Filter_Comparison.svg'));

%% ========================================================================
%  7. Generate background-subtracted fluorescence traces
%  ========================================================================
FbgSub_original = F - Fbg;
FbgSub_lowpass = F - Fbg_lowpass;
FbgSub_sgolay = F - Fbg_sgolay;
FbgSub_moving_avg = F - Fbg_moving_avg;

% Plotting the original F and the background-subtracted signals
figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.7]); % Wide figure

% Subplot 1: Original F
subplot(5,1,1);
plot(tvec, F, 'k', 'LineWidth', 1.5);
title('Original Fluorescence Signal (F)');
ylabel('Fluorescence (A.U.)');
grid on;

% Subplot 2: F - Original Background
subplot(5,1,2);
plot(tvec, FbgSub_original, 'm', 'LineWidth', 1.5);
title('F - Original Background');
ylabel('Fluorescence (A.U.)');
grid on;

% Subplot 3: F - Low-Pass Filtered Background
subplot(5,1,3);
plot(tvec, FbgSub_lowpass, 'r', 'LineWidth', 1.5);
title('F - Low-Pass Filtered Background');
ylabel('Fluorescence (A.U.)');
grid on;

% Subplot 4: F - Savitzky-Golay Filtered Background
subplot(5,1,4);
plot(tvec, FbgSub_sgolay, 'g', 'LineWidth', 1.5);
title('F - Savitzky-Golay Filtered Background');
ylabel('Fluorescence (A.U.)');
grid on;

% Subplot 5: F - Moving Average Filtered Background
subplot(5,1,5);
plot(tvec, FbgSub_moving_avg, 'b', 'LineWidth', 1.5);
title('F - Moving Average Filtered Background');
xlabel('Time (s)');
ylabel('Fluorescence (A.U.)');
grid on;

% Optimize the layout
sgtitle('Comparison of Background Subtracted Fluorescence Signals');

% Save the background subtracted figure
saveas(gcf, fullfile(pat, 'Background_Subtraction_Comparison.svg'));

%% ========================================================================
%  8. Save background-filter comparison variables
%  ========================================================================
save(fullfile(pat, 'Filtered_and_Subtracted_Fluorescence.mat'), ...
    'Fbg', 'Fbg_lowpass', 'Fbg_sgolay', 'Fbg_moving_avg', ...
    'FbgSub_original', 'FbgSub_lowpass', 'FbgSub_sgolay', 'FbgSub_moving_avg');

%% ========================================================================
%  9. Split field of view into a 4 x 4 ROI grid
%  ========================================================================

% Divide Idsm into 16 regions (4x4 grid)
rows = size(Idsm, 1);
cols = size(Idsm, 2);
nRegions = 4;
regionRows = floor(rows / nRegions);
regionCols = floor(cols / nRegions);

FbgSub_all = zeros(numFrames, nRegions, nRegions); % Pre-allocate FbgSub for all regions

for i = 1:nRegions
    for j = 1:nRegions
        regionMask = false(size(Idsm));
        regionMask((i-1)*regionRows+1:i*regionRows, (j-1)*regionCols+1:j*regionCols) = true;

        F_region = zeros(numFrames, 1); % Pre-allocate raw fluorescence trace for the region
        for k = 1:numFrames
            frame = Ids(:, :, k); % kth image frame
            F_region(k) = mean(frame(regionMask)); % Mean signal within the region mask of frame k
        end
        
        FbgSub_region = F_region; % Subtract background fluorescence
        FbgSub_all(:, i, j) = FbgSub_region; % Store FbgSub for the region
    end
end

regionLabels = cell(nRegions*nRegions, 1);

figure22Filename = fullfile(pat, 'ROI_heatmap.svg');
H22 = figure;

imagesc(Idsm); % Display the mean projection image
axis('off');
axis('image');
colormap(hot);
colorbar;
title('Mean Projection with Regions');

hold on;
for i = 1:nRegions
    for j = 1:nRegions
        x = (j-1) * regionCols + 1;
        y = (i-1) * regionRows + 1;
        rectangle('Position', [x, y, regionCols, regionRows], 'EdgeColor', 'w', 'LineWidth', 1.5);
        text(x + regionCols/2, y + regionRows/2, sprintf('(%d, %d)', i, j), 'Color', 'w', 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
end

saveas(H22, figure22Filename);

% Plot the FbgSub traces for each region
figure23Filename = fullfile(pat, 'FbgSub_traces_16regions.svg');
H23 = figure('position', [0, 0, 1600, 800]);

for i = 1:nRegions
    for j = 1:nRegions
        subplot(nRegions, nRegions, (i-1)*nRegions + j);
        plot(tvec, FbgSub_all(:, i, j), 'k');
        title(sprintf('Region (%d, %d)', i, j));
        if i == nRegions
            xlabel('Time (sec)');
        end
        if j == 1
            ylabel('F - F_{bg}');
        end
    end
end

saveas(H23, figure23Filename);

%% ========================================================================
%  10. Extract regional GRAB-ACh traces from the ROI grid
%  ========================================================================

% Settings
nRegions = 4;
[rows, cols] = size(Idsmbw);
regionRows = floor(rows / nRegions);
regionCols = floor(cols / nRegions);

FbgSub_all = nan(numFrames, nRegions, nRegions);  % Use NaN for excluded regions
dFF_all = nan(numFrames, nRegions, nRegions);     % Pre-allocate dF/F

for i = 1:nRegions
    for j = 1:nRegions
        % Define region bounds
        rowIdx = (i-1)*regionRows+1 : i*regionRows;
        colIdx = (j-1)*regionCols+1 : j*regionCols;

        % Mask: intersect region with GRAB-ACh mask
        regionMask = false(size(Idsmbw));
        regionMask(rowIdx, colIdx) = true;
        regionMask = regionMask & Idsmbw;  % Only include GRAB-ACh-positive pixels

        % Skip if regionMask is empty (e.g., no GRAB-ACh in this region)
        if ~any(regionMask(:))
            continue
        end

        % Extract raw F for region
        F_region = zeros(numFrames, 1);
        for k = 1:numFrames
            frame = Ids(:, :, k);
            F_region(k) = mean(frame(regionMask));
        end

        % Subtract global background
        FbgSub_region = F_region - Fbg;

        % Compute dF/F
        F0 = prctile(FbgSub_region, 5);
        dFF_region = (FbgSub_region - F0) / F0;

        % Store
        FbgSub_all(:, i, j) = FbgSub_region;
        dFF_all(:, i, j) = dFF_region;
    end
end

% Choose a representative image (e.g., mean projection of Ids)
meanImg = mean(Ids, 3);

% Display the image
figure;
imagesc(meanImg);
colormap(gray);
axis image off;
hold on;

% Overlay grid
[rows, cols] = size(Idsmbw);
regionRows = floor(rows / nRegions);
regionCols = floor(cols / nRegions);

for i = 1:nRegions
    for j = 1:nRegions
        % Define region bounds
        rowIdx = (i-1)*regionRows+1 : i*regionRows;
        colIdx = (j-1)*regionCols+1 : j*regionCols;

        % Only draw box if the region contains GRAB-ACh signal
        regionMask = false(size(Idsmbw));
        regionMask(rowIdx, colIdx) = true;
        regionMask = regionMask & Idsmbw;

        if any(regionMask(:))
            % Draw rectangle
            rectangle('Position', [colIdx(1), rowIdx(1), regionCols, regionRows], ...
                      'EdgeColor', 'r', 'LineWidth', 1.5);
            % Label it with its (i,j)
            text(colIdx(1)+2, rowIdx(1)+12, sprintf('(%d,%d)', i,j), ...
                 'Color', 'y', 'FontSize', 10, 'FontWeight', 'bold');
        else
            % Optional: outline in gray to show empty region
            rectangle('Position', [colIdx(1), rowIdx(1), regionCols, regionRows], ...
                      'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');
        end
    end
end

title('ROI grid overlaid on mean image');



%% ---- Display GRAB-ACh and background signal ROIs ---- % 
figure1Filename = fullfile(pat, 'heatmap.svg');
sz = get(0,'screensize');
H1 = figure('position',[1,sz(4)/2,sz(3),sz(4)/2]); 
movegui(H1, 'center')     

ax1 = subplot(1,3,1); 
imagesc(Idsm)                                       % display mean projection  
axis('off'), axis('image'), colorbar, colormap(ax1,hot) 
title('mean projection') 

ax2 = subplot(1,3,2); 
imagesc(Idsmbw),                                    % display binarised image 
axis('off'), axis('image'), colorbar, colormap(ax2,gray) 
title('GRAB-ACh area') 

ax3 = subplot(1,3,3); 
imagesc(IdsmBg),                                    % display binarised image 
axis('off'), axis('image'), colorbar, colormap(ax3,gray) 
title('background area') 

if exist('regions','var') && exist('roiX','var')    % if GRAB-ACh regions have been generated manually
    H2 = figure('color',[0.8 0.8 0.8]);                                    % new figure
    imagesc(Idsm)                                   % display mean projection
    axis('off'), axis('image'), colorbar, colormap(hot) 
    title('mean projection') 
    hold on
    plot(round(roiX),round(roiY),'m','linewidth',1.5)
    for j = 1:numel(regions)
        plot(regions{j}(:,2),regions{j}(:,1),'w','linewidth',1.5)
    end
    legend('Background area','GRAB-ACh area','location','southoutside','box','off')
end
saveas(H1, figure1Filename);
                                          % subtract background fluorescence
%%
if exist('Idsmbw','var') && exist('IdsmBg','var')
    H3 = figure('color',[0.8 0.8 0.8]); 
    imagesc(Idsm)                                   % display mean projection
    axis('off'), axis('image'), colorbar, colormap(hot)
    title('mean projection with mask outlines') 
    hold on

    % --- Overlay GRAB-ACh mask outline in white ---
    bwACh = bwboundaries(Idsmbw);
    for k = 1:length(bwACh)
        boundary = bwACh{k};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5)
    end

    % --- Overlay background mask outline in magenta ---
    bwBg = bwboundaries(IdsmBg);
    for k = 1:length(bwBg)
        boundary = bwBg{k};
        plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 1.5)
    end

    legend('GRAB-ACh area','Background area', ...
           'location','southoutside','box','off')
end

%%
Ffilt = FbgSub; % quick way to remove the high pass filter on the data for visualisation - comment out to keep high pass filter

F0 = median(Ffilt);                                 % baseline (F0) fluorescence
dFF = (Ffilt-F0)/F0;                                % delta F / F0 trace

% Plot dFF
figure;
plot(Ffilt, 'b', 'LineWidth', 1.5); % Plot dFF in blue
hold on;
yline(F0, 'r--', '5th Percentile F0', 'LineWidth', 1.5); % Indicate F0 with a red dashed line
xlabel('Time (frames)');
ylabel('dF/F');
title('dF/F Trace with 5th Percentile Baseline');
legend('dF/F', '5th Percentile F0');
grid on;
hold off;

%%
sz = get(0,'screensize');
% ---- display all fluorescence traces ---- %
figure2Filename = fullfile(pat, 'fluo_traces.svg');
H2 = figure('position',[1,1,sz(3)/2,sz(4)]); 
movegui(H2,'center')   

axA = subplot(5,1,1);
plot(tvec,F,'b')                                    % display raw fluorescence
ylabel('pix intensity')
set(gca,'xtick',[])
title('Raw fluorescence')

axB = subplot(5,1,2);
plot(tvec,Fbg,'r')                                  % display background fluorescence
ylabel('pix intensity')
set(gca,'xtick',[])
title('Background fluorescence')

axC = subplot(5,1,3);
plot(tvec,FbgSub,'g')                               % display background subtracted fluorescence
ylabel('pix intensity')
set(gca,'xtick',[])
title('Background subtracted fluorescence')

axD = subplot(5,1,4);
plot(tvec,Ffilt,'m')                                % display slow fluctuation filtered trace
ylabel('pix intensity')
set(gca,'xtick',[])
title('Slow filtered fluorescence')

axE = subplot(5,1,5);
plot(tvec,dFF,'k')                                  % display deltaF/F
ylabel('\DeltaF/F_0')
title('\DeltaF/F_0 trace')

xlabel('Time (sec)')
linkaxes([axA axB axC axD axE],'x')                 % link plot axes
set(gca,'xlim',[min(tvec) max(tvec)]);
saveas(H2, figure2Filename);

%%
sz = get(0,'screensize');

% ---- display fluorescence traces with background-subtracted and SG filtered dF/F ---- %
figure2Filename = fullfile(pat, 'fluo_traces_full.svg');
H2 = figure('position',[1,1,sz(3)/2,sz(4)]); 
movegui(H2,'center')   

% --- Raw fluorescence ---
axA = subplot(5,1,1);
plot(tvec,F,'b')                                    
ylabel('Pix intensity')
set(gca,'xtick',[])
title('Raw fluorescence')

% --- Background fluorescence ---
axB = subplot(5,1,2);
plot(tvec,Fbg,'r')                                  
ylabel('Pix intensity')
set(gca,'xtick',[])
title('Background fluorescence')

% --- Background-subtracted fluorescence ---
axC = subplot(5,1,3);
plot(tvec,FbgSub,'g')                              
ylabel('Pix intensity')
set(gca,'xtick',[])
title('Background-subtracted fluorescence')

% --- Delta F/F ---
axD = subplot(5,1,4);
plot(tvec,dFF,'k')                                  
ylabel('\DeltaF/F_0')
set(gca,'xtick',[])
title('Delta F/F')

% --- Savitzky-Golay filtered delta F/F ---
Fs = 1/mean(diff(tvec));        % sampling rate in Hz
win_samples = round(2 * Fs);    % 2-second window
if mod(win_samples,2)==0        % sgolayfilt requires odd window size
    win_samples = win_samples + 1;
end
dFF_sg = sgolayfilt(dFF, 3, win_samples);

axE = subplot(5,1,5);
plot(tvec,dFF_sg,'m')                               
ylabel('\DeltaF/F_0 (SG filtered)')
xlabel('Time (sec)')
title('Savitzky-Golay filtered Delta F/F')

% --- Link x-axes ---
linkaxes([axA axB axC axD axE],'x')                 
set(gca,'xlim',[min(tvec) max(tvec)]);

saveas(H2, figure2Filename);


%% ========================================================================
%  13. Load VR behaviour and align to imaging
%  ========================================================================
% -------------------------------------------------------------------------
% Load VR behaviour
% -------------------------------------------------------------------------

if vrWithimaging == 1
    fnVR = ['20' exptDate];
else
    if exptNo == 1
        fnVR = ['20' exptDate];
    else
        fnVR = ['20' exptDate '_' num2str(exptNo-1)];
    end
end

[vars,V,Y,t,licks,reward] = load_VR_data( ...
    behaviourRoot, ...
    animalID, ...
    fnVR, ...
    'Plot', false);

% generate indices denoting imaging start and end times
% Calculate the row indices for the selected session_
start_row = (imgNum - 1) * 2 + 1; % Calculate the start row index
end_row = (imgNum - 1) * 2 + 2;   % Calculate the end row index

% Get the start and end times for the selected session from vars.times.laser_times
session_start_time = vars.times.laser_times(start_row, 2);
session_end_time = vars.times.laser_times(end_row, 2);
session_indices = find(t >= session_start_time & t <= session_end_time);
t_session = t(session_indices);

% create VR behaviour variables within the laser_times imaging window
tNorm = t_session-min(t_session);                                   % normalise t so the time series starts at 0
% tNorm = tNorm(session_indices);
V_session = V(session_indices);
licks_session = licks(session_indices);
reward_session = reward(session_indices);
Y_session = Y(session_indices);

% create a downsampled velocity vector that matches the imaging sampling frequency
Vds = nan(length(tvec),1);                          % pre-allocate downsampled velocity array
for k = 1:length(tvec)
    [~,ind] =  min(abs(tNorm-tvec(k)));             % find the timestamp in tNorm that most closely matches the kth imaging timestamp (i.e. in tvec) 
    Vds(k) = V_session(ind);                                % add the velocity for the closest matched timestamp to the downsampled velocity array   
end
VdsSm = smooth(Vds,Fsds);                           % smooth the velecity with a 1 second moving average filter 
% VdsSm2 = smooth(Vds, Fsds*0.5);
%%% create a downsampled licking vector matching imaging sampling frequency
% Define the width of the time window for averaging licks (in seconds)
window_width = 1.0; % You can adjust this value as needed
% Calculate the number of data points in each window based on the window width and sampling frequency
window_size = round(window_width * Fsds);
% Pre-allocate the downsampled lick array
licks_ds = nan(length(tvec), 1);

for k = 1:length(tvec)
    [~, ind] = min(abs(tNorm - tvec(k))); % Find the timestamp in tNorm that most closely matches the kth imaging timestamp (i.e., in tvec)
    
    % Define the start and end indices for the time window
    start_idx = max(1, ind - window_size + 1);
    end_idx = min(length(licks_session), ind);
    
    % Calculate the average number of licks within the window
    licks_ds(k) = sum(licks_session(start_idx:end_idx)) / (end_idx - start_idx + 1);
end


% create a downsampled y (position) vector
Y_ds = nan(length(tvec), 1);
for k = 1:length(tvec)
    [~, ind] = min(abs(tNorm - tvec(k)));
    Y_ds(k) = Y_session(ind);
end


%%
% Set the threshold for lap transitions
threshold = 100;

% Identify indices where Y_session decreases significantly (lap transitions)
lap_transitions = find(diff(Y_session) < -threshold);
lap_transitions = [1; lap_transitions; length(Y_session)];  % Include the first and last index as lap boundaries

% Initialize downsampled reward vector
reward_ds = zeros(length(tvec), 1);

% Iterate through each lap
for lap = 1:length(lap_transitions)-1
    lap_start = lap_transitions(lap);
    lap_end = lap_transitions(lap+1);
    
    % Find the reward position in the current lap (relative index)
    reward_in_lap_relative = find(reward_session(lap_start:lap_end) == 1);
    
    if ~isempty(reward_in_lap_relative)
        % Find the first reward in this lap
        reward_in_lap_absolute = lap_start + reward_in_lap_relative(1) - 1;
        
        % Get the corresponding Y_session value for the reward in this lap
        reward_position = Y_session(reward_in_lap_absolute);
        
        % Restrict Y_ds to the current lap range for better precision
        lap_ds_start = find(tvec >= tNorm(lap_start), 1, 'first');
        lap_ds_end = find(tvec <= tNorm(lap_end), 1, 'last');
        
        if ~isempty(lap_ds_start) && ~isempty(lap_ds_end)
            % Find the closest Y_ds value to this reward_position within the current lap's Y_ds range
            [~, closest_row] = min(abs(Y_ds(lap_ds_start:lap_ds_end) - reward_position));
            closest_row = lap_ds_start + closest_row - 1;  % Convert to the global index in Y_ds
            
            % Set the reward_ds value to 1 at the closest row, only for this lap
            reward_ds(closest_row) = 1;
        end
    end
end

% Plot Y_session with reward_session overlaid
figure;
subplot(2, 1, 1);  % First plot: Original Y_session with reward_session
plot(tNorm, Y_session, 'b', 'LineWidth', 1.5);  % Plot Y_session
hold on;
reward_indices = find(reward_session == 1);
plot(tNorm(reward_indices), Y_session(reward_indices), 'ro', 'MarkerFaceColor', 'r');  % Overlay reward locations
title('Original Y\_session with Reward Locations');
xlabel('Time');
ylabel('Y\_session');
legend('Y\_session', 'Reward Locations');
grid on;
hold off;

% Plot Y_ds with reward_ds overlaid
subplot(2, 1, 2);  % Second plot: Downsampled Y_ds with reward_ds
plot(tvec, Y_ds, 'b', 'LineWidth', 1.5);  % Plot Y_ds
hold on;
reward_ds_indices = find(reward_ds == 1);
plot(tvec(reward_ds_indices), Y_ds(reward_ds_indices), 'ro', 'MarkerFaceColor', 'r');  % Overlay reward locations
title('Downsampled Y\_ds with Reward Locations');
xlabel('Time');
ylabel('Y\_ds');
legend('Y\_ds', 'Reward Locations');
grid on;
hold off;


%% ========================================================================
%  14. Plot aligned GRAB-ACh and behaviour traces
%  ========================================================================
figure3Filename = fullfile(pat, 'GRAB_VR.svg'); % You can choose a suitable file extension, e.g., .png
H3 = figure('Position', [100, 100, 800, 600]);

% Create the first subplot for tvec and F
ax4 = subplot(5,1,1);
% plot(tvec, F, 'g')                                    % Display F trace
ylabel('F')
set(gca, 'xtick', [])
title('Fluorescence, \DeltaF/F_0, Velocity, Licking, and Reward Traces ')  % Add a title if needed
xlim(ax4, [min(tvec), max(tvec)]);

% Create the second subplot for dFF
ax1 = subplot(5,1,2);
plot(tvec, dFF, 'k')                                  % Display GRAB-ACh trace
ylabel('\DeltaF/F_0')
set(gca, 'xtick', [])
xlim(ax1, [min(tvec), max(tvec)]);

% Create the third subplot for VdsSm (running speed)
ax2 = subplot(5,1,3);
plot(tvec, VdsSm, 'r')                                % Display velocity trace
ylabel({'Running speed'; '(cm/sec)'})
set(gca, 'xtick', [])
xlim(ax2, [min(tvec), max(tvec)]);

% Create the fourth subplot for Licks_ds (downsampled licks)
ax3 = subplot(5,1,4);
plot(tvec, licks_ds, 'b')                            % Display downsampled lick data
ylabel('Licks')
ylim([0 2])
set(gca, 'xtick', [])
xlim(ax3, [min(tvec), max(tvec)]);

% Create the fifth subplot for reward_ds
ax5 = subplot(5,1,5);
plot(tvec, reward_ds, 'm')                            % Display reward_ds data
ylabel('Reward')
xlabel('Time (sec)')
xlim(ax5, [min(tvec), max(tvec)]);

% Link all subplots along the x-axis
linkaxes([ax1, ax2, ax3, ax4, ax5], 'x')
saveas(H3, figure3Filename);

% ---- Compare slow Load behaviour data ---- %
indRest = VdsSm < 1;                                % indices of epochs where mouse is resting (speed < 1 cm/sec)
dFFRest = mean(dFF(indRest));                       % GRAB-ACh fluorescence when mouse is resting

indRun = VdsSm > 5;                                 % indices of epochs where mouse is running (speed > 5 cm/sec)
dFFRun = mean(dFF(indRun));                         % GRAB-ACh fluorescence when mouse is running

% plot rest vs. run comparison
% traces 
figure4Filename = fullfile(pat, 'rest_vs_run.svg'); % You can choose a suitable file extension, e.g., .png
H4 = figure; 
plot([1 2],[dFFRest dFFRun],'bo')
hold on
plot([1 2],[dFFRest dFFRun],'b')
set(gca,'xlim',[0.5 2.5],'xtick',[1 2],'xticklabel',{'Rest','Run'})
ylabel('GRAB-ACh \DeltaF/F_0')
saveas(H4, figure4Filename);

% ---- Perform spectral analysis using continuous wavelent transform ---- %
nVoices = 30;                                       % voices per octave - frequency resolution
[wt,freq] = cwt(dFF,'amor',Fsds,'VoicesPerOctave',nVoices); % wavelet transform (Morlet wavelet function)
wta = abs(wt);                                      % wavelket transform amplitude  

% display data
figure5Filename = fullfile(pat, 'freq.svg'); % You can choose a suitable file extension, e.g., .png
H5 = figure;
ax1 = subplot(4,1,1);
restMarker = nan(length(tvec),1);
restMarker(indRest) = -2;
plot(tvec,restMarker,'r','linewidth',2)             % display rest epochs
hold on
runMarker = nan(length(tvec),1);
runMarker(indRun) = -4;                             
plot(tvec,runMarker,'g','linewidth',2)              % display run epochs
L = legend('Rest','Run','location','northeast','box','off');
set(L,'AutoUpdate','off')
plot(tvec,VdsSm,'b')                                % display running velocity
set(gca,'xtick',[])
ylabel('cm/sec')                               

ax2 = subplot(4,1,2);
plot(tvec,dFF,'k')                                  % display deltaF/F
set(gca,'xtick',[])
ylabel('\DeltaF/F_0')

ax3 = subplot(4,1,3:4);
imagesc(tvec,freq,wta)                              % display wavelet transform
set(gca,'ydir','normal')
set(gca,'yscale','log')
set(gca,'ylim',[0.05 15])
set(gca,'yminortick','on')
ylabel('Frequency (Hz)')
colormap(hot)
ax3.YColor = 'red';
ax3.TickLength = [0.03 0.03];

xlabel('Time (sec)')
linkaxes([ax1 ax2 ax3],'x')                         % link plot axes
set(gca,'xlim',[min(tvec) max(tvec)])
saveas(H5, figure5Filename);

%% ========================================================================
%  15. Save processed GRAB-ACh and behaviour variables
%  ========================================================================

% Create the save file name
saveFileName = 'data.mat';

% Determine the save file path using the 'pat' variable
saveFilePath = fullfile(pat, saveFileName);

% Specify the variables you want to save
save(saveFilePath, 'F', 'FbgSub', 'dFF', 'tvec', 'Fsds', 'VdsSm', 'licks_ds', 'reward_ds', 'Y_ds', 'Fs', 'numFrames', 'pat', 'dFFRun', 'dFFRest', 'tNorm');

% with no behaviour
% save(saveFilePath, 'F', 'FbgSub', 'dFF', 'tvec', 'Fsds', 'Fs', 'numFrames', 'pat');
