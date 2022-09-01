% Load up a fully sampled 8-channel K-space example
% kspace is full 2D data, [nRO, nPE, nCh]
load brain_data_8ch_noisy.mat
[nRO,nPE,nCh] = size(kspace);

% Simulate undersampling by only keeping Rth line
R = 8;
kspUndersamp = zeros(nRO,nPE,nCh);
kspUndersamp(:,1:R:end,:) = kspace(:,1:R:end,:);

% Display kspace for channel 1
figure(1)
subplot(1,2,1)
imagesc(log(abs(kspace(:,:,1))));
title('log(abs(ksp)) for Ch #1, full')

subplot(1,2,2)
imagesc(log(abs(kspUndersamp(:,:,1))));
title(sprintf('log(abs(ksp)) for Ch #1, undersampled %dx', R))
linkaxes 
zoom on

%% Sensitivity profiles
% Convert to image space. img_mc = multi-channel complex image
img_mc = fftshift(ifft2(ifftshift(kspace))); 

% Calculate sensitivity profiles using root sum-of-squared denominator
% with fully-sampled data. 
img_RSOS = sqrt(sum(img_mc .* conj(img_mc), 3));
sensMap = zeros(nRO,nPE,nCh);
for iCoil=1:nCh
    sensMap(:,:,iCoil) = img_mc(:,:,iCoil)./img_RSOS;
end

% Display sensitivity profiles
figure(2)
for iCoil = 1:nCh
    subplot(2, nCh, iCoil);
    imagesc(abs(sensMap(:,:,iCoil)))
    title(sprintf('abs(%d)', iCoil));

    subplot(2, nCh, iCoil+nCh);
    imagesc(angle(sensMap(:,:,iCoil)))
    title(sprintf('angle(%d)', iCoil));
end

%% SENSE
% Run the SENSE reconstruction on the undersampled data
imgRecon = cgSENSE(sensMap, kspUndersamp);

% Evaluate difference relative to our reference RSOS scan
diff = img_RSOS - imgRecon;
mae = mean(abs(diff(:)));
rmse = sqrt(mean(diff(:).^2)); 
fprintf('mean absolute error %f, rmse %f\n', mae, rmse);

% Display the images
figure(3)
colormap parula
clim = [0 3]; % limits of intensity to display
subplot(2,2,1)
imagesc(img_RSOS, clim)
title('RSOS recon of fully sampled data');

subplot(2,2,2)
imagesc(imgRecon, clim)
title(sprintf('Reconstructed image, R=%d', R))

subplot(2,2,3)
imagesc(abs(diff)*100, clim)
title('abs(difference) x100')

% This is so you can zoom in the same across all images
linkaxes; zoom on




%% PJB hack
figure(12)
colormap gray
img = imgaussfilt(imgRecon, 1);
img = imresize(img, 2, 'bicubic');    
imagesc(img, [.7 1.8]);
zoom on


