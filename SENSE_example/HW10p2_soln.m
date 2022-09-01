% Load up a fully sampled 8-channel K-space example
% kspace is full 2D data, [nRO, nPE, nCh]
load brain_data_8ch_noisy.mat

% Convert to image space
img_mc = fftshift(ifft2(ifftshift(kspace)));
[nRO, nPE, nCh] = size(kspace);
R = 2;

% Calculate sensitivity maps.
sens_8ch = calculate_sensitivity_maps(img_mc);

% For my gold-standard reference, I'll use the RSOS of all 8 channels,
% fully sampled.
imgRSOS_8ch = sqrt(sum(abs(fftshift(ifft2(ifftshift(kspace)))).^2,3));
fprintf('\n\n')

% Loop over all 2nd channels
mae = zeros(1,8);
for iCoil = 2:8
    
    ch_select = [1, iCoil];
    
    % Simulate undersampling by only keeping Rth line, and only the 2
    % channels:
    kspUndersamp = zeros(nRO, nPE, 2);
    kspUndersamp(:,1:R:end,:) = kspace(:,1:R:end, ch_select);
    
        
    % Reconstruct
    imgRecon_2ch = cgSENSE(sens_8ch(:,:,ch_select), kspUndersamp);
    % Alterante interpretation: see solution
    %sens_2ch = calculate_sensitivity_maps(img_mc(:,:, ch_select));
    %imgRecon_2ch = cgSENSE(sens_2ch, kspUndersamp); 
  
    
    diff = imgRSOS_8ch - imgRecon_2ch;
    % Alterante interpretation: see solution    
    %imgRSOS_2ch = sqrt(sum(abs(fftshift(ifft2(ifftshift(kspace(:,:,ch_select))))).^2,3));
    %diff = imgRSOS_2ch - imgRecon_2ch;

    % Measure error, save the value
    mae(iCoil) = mean(abs(diff(:)));
    
    fprintf('Ch %d: mean absolute error %f\n', iCoil, mae(iCoil));
    
end

%
figure(1)
bar(2:8, mae(2:8))
xlabel('2nd channel')
ylabel('MAE')

%%

% Evaluate reconstructions with optimized display
figure(12)
colormap gray
subplot(2,4,1);
plot_image_optimized(imgRSOS_8ch)
title('Reference')

for iCoil = 2:8
    ch_select = [1, iCoil];
    kspUndersamp = zeros(nRO, nPE, 2);
    kspUndersamp(:,1:R:end,:) = kspace(:,1:R:end, ch_select);
    imgRecon_2ch = cgSENSE(sens_8ch(:,:,ch_select), kspUndersamp);
    
    subplot(2,4,iCoil);
    plot_image_optimized(imgRecon_2ch);
    title(sprintf('Ch 1,%d', iCoil));
end
linkaxes

%%

function plot_image_optimized(img)
img = imgaussfilt(img, 1);
img = imresize(img, 2, 'bicubic');
imagesc(img(100:250, 200:350), [0 1.5]);
%imagesc(img, [0 1.5]);
axis equal; axis tight
zoom on
end


% Calculate sensitivity map using RSOS for denominator
function sens = calculate_sensitivity_maps(img_mc)

[nRO, nPE, nCh] = size(img_mc);
img_RSOS = sqrt(sum(abs(img_mc).^2,3));

sens = zeros(nRO,nPE,nCh);
for iCoil=1:nCh
    sens(:,:,iCoil) = img_mc(:,:,iCoil)./img_RSOS;
end

end



