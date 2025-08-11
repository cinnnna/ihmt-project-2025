%% Step 1: Load and calculate ihMT percentage
M0 = double(dicomread('M0_185_1.dcm'));
ihMT = double(dicomread('ihMT_188_2.dcm'));

% Apply median filter to ihMT for denoising
ihMT_denoised = medfilt2(ihMT, [3 3]);

% Calculate ihMT ratio as percentage
ihMTR_percent = (ihMT_denoised ./ M0) * 100;

% Handle any inf/nan values (where M0 might be zero)
ihMTR_percent(isinf(ihMTR_percent)) = 0;
ihMTR_percent(isnan(ihMTR_percent)) = 0;

%% Step 2: Determine appropriate window for percentage data
% Check the range of percentage values
brain_mask = M0 > prctile(M0(:), 10); % Simple brain mask
brain_values = ihMTR_percent(brain_mask);

% Set window for percentage - typical ihMT is 0-10%
window_percent = [0, 10]; % Adjust based on your data
% Or use percentile-based:
% window_percent = [prctile(brain_values, 1), prctile(brain_values, 99)];

%% Step 3: Apply windowing and normalize
ihMTR_windowed = ihMTR_percent;
ihMTR_windowed(ihMTR_windowed < window_percent(1)) = window_percent(1);
ihMTR_windowed(ihMTR_windowed > window_percent(2)) = window_percent(2);

% Normalize to 0-1 for display
ihMTR_norm = mat2gray(ihMTR_windowed, window_percent);

%% Step 4: Create the blue (jet) colored image
figure('Color', 'white', 'Position', [100 100 800 600]);

% Display image
imshow(ihMTR_norm);
colormap(jet); % Blue colormap as requested

% Add professional colorbar with percentage units
c = colorbar;
c.Label.String = 'ihMT Ratio (%)';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';
c.FontSize = 12;

% Update colorbar ticks to show percentage values
num_ticks = 6;
tick_positions = linspace(0, 1, num_ticks);
c.Ticks = tick_positions;
tick_values = linspace(window_percent(1), window_percent(2), num_ticks);
c.TickLabels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);

% Add title
title('ihMT Ratio Map', 'FontSize', 18, 'FontWeight', 'bold');
axis off;

% Save high-resolution image
exportgraphics(gcf, 'ihMTR_percentage_jet.png', 'Resolution', 300);

%% Optional: Create comparison figure
figure('Color', 'white', 'Position', [50 50 1600 600]);

% Original ihMT in arbitrary units
subplot(1,3,1);
imshow(mat2gray(ihMT_denoised, [300, 7500]));
colormap(gca, jet);
c1 = colorbar;
c1.Label.String = 'ihMT ratio M_z/M_0 (%)';
title('ihMT ratio (%)', 'FontSize', 14);

% M0 reference
subplot(1,3,2);
imshow(M0, []);
colormap(gca, gray);
c2 = colorbar;
c2.Label.String = 'M0 Signal (a.u.)';
title('M0 Reference', 'FontSize', 14);

% ihMT percentage
subplot(1,3,3);
imshow(ihMTR_norm);
colormap(gca, jet);
c3 = colorbar;
c3.Label.String = 'ihMT Ratio (%)';
% Update ticks for percentage
c3.Ticks = tick_positions;
c3.TickLabels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);
title('ihMT Ratio (%)', 'FontSize', 14);

sgtitle('ihMT Processing: From Arbitrary Units to Percentage', 'FontSize', 16);