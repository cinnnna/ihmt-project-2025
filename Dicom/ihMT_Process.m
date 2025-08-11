% Load your image and check the distribution
img = double(dicomread('ihMT_188_2.dcm'));

% Get brain tissue statistics (exclude background)
brain_mask = img > 100;  % Rough brain mask
brain_pixels = img(brain_mask);

% Calculate percentiles
percentiles = prctile(brain_pixels, [1, 5, 25, 50, 75, 95, 99]);
fprintf('Brain tissue percentiles:\n');
fprintf('1%%: %.0f\n', percentiles(1));
fprintf('5%%: %.0f\n', percentiles(2));
fprintf('25%%: %.0f\n', percentiles(3));
fprintf('50%% (median): %.0f\n', percentiles(4));
fprintf('75%%: %.0f\n', percentiles(5));
fprintf('95%%: %.0f\n', percentiles(6));
fprintf('99%%: %.0f\n', percentiles(7));

% Your current window coverage
your_window = [300, 7500];
coverage = sum(brain_pixels >= your_window(1) & brain_pixels <= your_window(2)) / length(brain_pixels) * 100;
fprintf('\nYour window [%d, %d] covers %.1f%% of brain tissue\n', your_window(1), your_window(2), coverage);


%% Step-by-Step Processing for ihMT Image

% STEP 1: Load the original image
img_original = double(dicomread('ihMT_188_2.dcm'));
fprintf('Step 1: Image loaded. Size: %dx%d\n', size(img_original));

% STEP 2: Apply median filter to remove spike noise
img_denoised = medfilt2(img_original, [3 3]);
fprintf('Step 2: Median filter applied\n');

% STEP 3: Apply your window settings
img_windowed = img_denoised;  % Work on the denoised image
img_windowed(img_windowed < 300) = 300;  % Set minimum
img_windowed(img_windowed > 7500) = 7500;  % Set maximum
fprintf('Step 3: Window levels applied [1112, 5556]\n');

% STEP 4: Normalize to 0-1 range for display
img_normalized = mat2gray(img_windowed, [300, 5034]);
fprintf('Step 4: Normalized to [0,1] range\n');

% STEP 5: Optional - Apply mild smoothing
img_final = imgaussfilt(img_normalized, 1);
fprintf('Step 5: Gaussian smoothing applied\n');

% STEP 6: Display the progression
figure('Position', [50 50 1800 400], 'Name', 'Processing Steps');

subplot(1,5,1);
imshow(img_original, []); 
title('1. Original', 'FontSize', 12);

subplot(1,5,2);
imshow(img_denoised, []); 
title('2. After Median Filter', 'FontSize', 12);

subplot(1,5,3);
imshow(img_windowed, [300, 7500]); 
title('3. After Windowing', 'FontSize', 12);

subplot(1,5,4);
imshow(img_normalized); 
title('4. After Normalization', 'FontSize', 12);

subplot(1,5,5);
imshow(img_final); 
title('5. Final (with smoothing)', 'FontSize', 12);

% STEP 7: Save the final result
imwrite(img_final, 'ihMT_processed_final.png');
fprintf('Step 7: Final image saved\n');