% Load up a fully sampled 8-channel K-space example
% kspace is full 2D data, [nRO, nPE, nCh]
load brain_data_8ch_noisy.mat
[nRO,nPE,nCh] = size(kspace);

% Simulate undersampling by only keeping Rth line
R = 2;
kspUndersamp = zeros(nRO,nPE,nCh);
kspUndersamp(:,1:R:end,:) = kspace(:,1:R:end,:);

% Make logical masks of which points were sampled
maskSkipped = logical(kspUndersamp == 0);
maskAcquired = ~maskSkipped;

% Calculate sensitivity profiles using root sum-of-squared denominator
% with fully-sampled data. If you are doing an accelerated acqusition this
% information would have to come from a separate scan
imgRSOS = sqrt(sum(abs(fftshift(ifft2(ifftshift(kspace)))).^2,3));
channelImages = fftshift(ifft2(ifftshift(kspace)));
sensMap = zeros(nRO,nPE,nCh);
for iCoil=1:nCh
    sensMap(:,:,iCoil) = channelImages(:,:,iCoil)./imgRSOS;
end

% Reconstruct with a subset of channels
ch_select = [1, 8];
nCh = size(ch_select, 2);
kspUndersamp = kspUndersamp(:, :, ch_select);
sensMap = sensMap(:,:,ch_select);

% Run the SENSE reconstruction on the undersampled data
imgRecon = cgSENSE(sensMap, kspUndersamp);

% Evaluate difference with our reference RSOS scan
diff = imgRSOS - imgRecon;
mae = mean(abs(diff(:)));
rmse = sqrt(mean(diff(:).^2)); 
fprintf('CH %d: mean absolute error %f, rmse %f\n', ch_select(2), mae, rmse);

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


% Display sensitivity profiles
figure(2)
for cdx = 1:nCh
    subplot(2, nCh, cdx);
    imagesc(abs(sensMap(:,:,cdx)))
    title(sprintf('abs(%d)', cdx));

    subplot(2, nCh, cdx+nCh);
    imagesc(angle(sensMap(:,:,cdx)))
    title(sprintf('angle(%d)', cdx));
end


% Display the images
figure(3)
clim = [0 5];
subplot(2,2,1)
imagesc(imgRSOS, clim)

title('RSOS recon of fully sampled data');

subplot(2,2,2)
imagesc(imgRecon, clim)
title(sprintf('Reconstructed image, R=%d', R))

subplot(2,2,3)
imagesc(abs(diff)*100, clim)
title('abs(difference) x100')
linkaxes 
zoom on



