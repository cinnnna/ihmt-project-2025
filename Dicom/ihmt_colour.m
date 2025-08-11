%% Load ihMT data
ihMT = double(dicomread('ihMT_188_2.dcm'));

% Apply median filter for denoising
ihMT_denoised = medfilt2(ihMT, [3 3]);

% Apply your window [300, 7500]
window = [300, 7500];
ihMT_windowed = ihMT_denoised;
ihMT_windowed(ihMT_windowed < window(1)) = window(1);
ihMT_windowed(ihMT_windowed > window(2)) = window(2);

% Normalize to 0-1 for display
ihMT_norm = mat2gray(ihMT_windowed, window);

%% Create and save the blue (jet) colored image
figure('Color', 'white', 'Position', [100 100 800 600]);
imshow(ihMT_norm);
colormap(jet);

% Add colorbar with percentage labels (0-100)
c = colorbar;
c.Label.String = 'ihMT Signal (%)';
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';
c.FontSize = 12;

% Change colorbar ticks to show 0-100 instead of 0-1
c.TickLabels = {'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'};

% Add title
%title('ihMT Contrast Map', 'FontSize', 18, 'FontWeight', 'bold');
%axis off;

% Save high-resolution image
exportgraphics(gcf, 'ihMT_jet_normalized.png', 'Resolution', 300);