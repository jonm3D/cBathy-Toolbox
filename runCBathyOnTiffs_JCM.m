%% Prep and Running CBathy
% This script processes a set of rectified image frames to estimate bathymetry (underwater depth) using the CBathy method.

% Clear previous data and figures (if plotting is enabled)
clear all;
close all;

%% Control Plotting

% Plotting Flag (Set to true to enable plotting, false to disable)
enablePlotting = false; % JCM plots

%% List Available Directories and Select One

% Base directory containing all sXXX... directories
baseDir = '/Users/jonathan/Desktop/frf_collects/';

% Get a list of directories matching the sXXX pattern
dirInfo = dir(fullfile(baseDir, 's*'));
dirNames = {dirInfo([dirInfo.isdir]).name};

% Display available directories and let the user select one
fprintf('Available directories:\n');
for i = 1:length(dirNames)
    fprintf('%d: %s\n', i, dirNames{i});
end
fprintf('%d: Process all directories\n', length(dirNames) + 1);

selectedIdx = input('Select a directory by entering its number: ');

% Validate the selection
if selectedIdx < 1 || selectedIdx > length(dirNames) + 1
    error('Invalid selection. Please restart the script and select a valid directory.');
end

if selectedIdx == length(dirNames) + 1
    % Process all directories
    for i = 1:length(dirNames)
        processDirectory(fullfile(baseDir, dirNames{i}), enablePlotting);
    end
else
    % Process the selected directory
    selectedDir = fullfile(baseDir, dirNames{selectedIdx});
    processDirectory(selectedDir, enablePlotting);
end

%% Function to process a single directory
function processDirectory(selectedDir, enablePlotting)
    fprintf('Processing directory: %s\n', selectedDir);
    framesDir = fullfile(selectedDir, 'frames');
    cd(framesDir);

    %% Load Image Frames and Define Grid

    % Get collection identifier
    id_pattern = 's\d{3}_\d{8}';

    % Use regular expression to find the match
    tokens = regexp(selectedDir, id_pattern, 'match');

    % Check if a match is found
    if ~isempty(tokens)
        collectId = tokens{1}; % Extract the whole string starting with sXXX...
    else
        collectId = ''; % Return an empty string if no match is found
    end

    % Get current date and time
    dateTimeStr = datestr(now, 'yyyymmdd_HHMMSS');

    % Load image file names
    dc = dir('*.tiff');

    % Extract the base name of the directory and parent directory
    [~, baseName, ~] = fileparts(fileparts(selectedDir));
    [parentDir, parentFolder, ~] = fileparts(fileparts(selectedDir));

    % Define the output directory with date and time suffix
    outputDir = fullfile(parentDir, baseName, ['output_' dateTimeStr]);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Define the full image coordinate grid (X, Y)
    % Images should already be registered/stabilized

    Xv = -700:1:1299;
    Yv = -300:1:1699;
    [X, Y] = meshgrid(Xv, Yv);

    % Define time resolution and time vector
    tRes = 1; % seconds or frames
    tTot = length(dc) / 30; % seconds
    t = 1/30:(1/30):tTot;

    % Define the region of interest (ROI) on the image grid
    res = 2;
    sx = 700;
    ex = 1800;
    sy = 300;
    ey = 1300;

    % Create meshgrid for the ROI
    [Xg, Yg] = meshgrid(Xv(sx:res:ex), Yv(sy:res:ey));

    %% Load and Preprocess Images

    % Load the first image (and flip vertically if needed) for sizing
    I = flipud(imread(fullfile(framesDir, dc(1).name)));
    Ig = I(sy:res:ey, sx:res:ex);

    % Preallocate time stack for efficiency
    TS = zeros([size(Ig), length(1:tRes:length(dc))]);

    % Ensure Parallel Computing Toolbox is available
    if isempty(gcp('nocreate'))
        parpool;
    end

    % Loop through all frames with the chosen time resolution using parfor
    parfor i = 1:tRes:length(dc)
        image = flipud(imread(fullfile(framesDir, dc(i).name)));
        TS(:, :, i) = adapthisteq(image(sy:res:ey, sx:res:ex));  % Apply adaptive histogram equalization
    end

    %% Visualize the Region of Interest (Conditional Plotting)

    if enablePlotting
        figure;
        pcolor(Xg, Yg, TS(:, :, 1));
        shading flat;
        xlabel('X [m]');
        ylabel('Y [m]');
        colormap(gray);
        hold on;
        rectangle('Position', [min(min(Xg)) min(min(Yg)) max(max(Xg))-min(min(Xg)) max(max(Yg))-min(min(Yg))], 'edgecolor', 'g');

        % Save the figure
        saveas(gcf, fullfile(outputDir, sprintf('RegionOfInterest_%s.fig', collectId)));
        saveas(gcf, fullfile(outputDir, sprintf('RegionOfInterest_%s.png', collectId)));
    end

    %% Prepare Data for CBathy

    tepoch_cut = t; %(t-datenum(1970,1,1)).*24*3600; %cbathy likes epoch time

    % Reshape Image Data for CBathy
    [r, c] = size(TS(:,:,1));
    [xindgrid, yindgrid] = meshgrid(1:c, 1:r);
    rowIND = yindgrid(:);
    colIND = xindgrid(:);

    counter = 0;
    data = [];  
    for i = 1:length(rowIND(:))
        counter = counter + 1;
        data(:, counter) = reshape(TS(rowIND(i), colIND(i), :), length(tepoch_cut), 1);
    end

    %% Set Up CBATHY

    xyz = [Xg(:) Yg(:) Xg(:) * 0];  
    argus02b; 
    bathy.params = params;  
    bathy.params.xyMinMax = [0, 900, 0, 1000]; 
    cam = ones(size(xyz, 1), 1);

    %% Run CBathy and Analyze Results

    % (Optional) Filter specific frequency bands for the analysis
    % bathy.params.fB = bathy.params.fB(4:10); 

    [f, G, bathy] = prepBathyInputShort(xyz, tepoch_cut, double(data), bathy);
    fprintf('Done preparing input for analysis, running analyzeBathyCollect...')
    bathy = analyzeBathyCollectShort(xyz, tepoch_cut, double(data), cam, bathy); 
    fprintf('Done analyzing data, generating plots...')

    %% Plot Results (Conditional Plotting)

    if enablePlotting
        % Plot the original image and the estimated bathymetry side-by-side
        f1 = figure;
        a = gca; 
        pcolor(Xg, Yg, TS(:, :, 1));
        colormap(a, "gray");
        shading flat;
        set(gca, 'ydir', 'normal');
        hold on;
        a2 = axes;
        linkaxes([a a2], 'xy');
        a2.Color = 'none';
        hold on;
        pcolor(bathy.xm, bathy.ym, -bathy.fCombined.h);
        shading flat;
        h = colorbar('south');
        caxis([-10 0]);
        h.Color = 'w';
        h.Label.String = 'Depth [m]';
        xlabel('X [m]');
        ylabel('Y [m]');

        % Save the figure
        saveas(f1, fullfile(outputDir, sprintf('Bathy_%s.fig', collectId)));
        saveas(f1, fullfile(outputDir, sprintf('Bathy_%s.png', collectId)));
        
        figure(13);
        i = 4;  
        ind = find(abs(f - params.fB(i)) == min(abs(f - params.fB(i))));
        h = scatter(xyz(:, 1), xyz(:, 2), 3, angle(G(ind, :)), 'filled'); 
        xlabel('X [m]');
        ylabel('Y [m]');

        % Save the figure
        saveas(gcf, fullfile(outputDir, sprintf('PhaseMap_%s.fig', collectId)));
        saveas(gcf, fullfile(outputDir, sprintf('PhaseMap_%s.png', collectId)));
    end

    %% Save Results (Optional)
    outputFileName = sprintf('cBathy2Hz%.0fs_%s.mat', tTot, collectId);
    save(fullfile(outputDir, outputFileName), "bathy", "data", "params", "f", "G", "Xg", "Yg", "TS", "xyz", "cam", "tepoch_cut", '-v7.3');
end
