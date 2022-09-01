%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data in format nRO, nPE, nCh
% This is an implementation of 2D conjugate gradient SENSE
% Written by Steen Moeller with adaptations by Jessica McKay and Pat Bolan
% Inputs:
%   sensMap:        complex sensitivity map [nRO, nPE, nCH]
%   kspUndersamp:   the undersampled kspace data. The matrix should have 
%                   the size of fully sampled data, [nRO, nPE, nCH], with
%                   zeros for unsampled values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [imgRecon, err] = cgSENSE(sensMap, kspUndersamp)

[nRO, nPE, nCh] = size(kspUndersamp);

maskSkipped = logical(kspUndersamp == 0);
maskAcquired = ~maskSkipped;
y_measured = kspUndersamp(maskAcquired);

E = @(x) sense_op(sensMap, x, nRO, nPE, nCh, maskAcquired);
ET = @(x) sense_op_t(sensMap, x, maskAcquired, nRO, nPE, nCh);
ATA = @(x) ET(E(x));
ATb = ET(y_measured);


inp = sum(ifft2(squeeze(kspUndersamp(:,:,:))),3);

[imm, err] = conjgrad(20, ATA, ATb, inp(:),kspUndersamp,maskAcquired);
imgRecon = abs(reshape(imm,[nRO nPE]));

return;


function [x, err] = conjgrad(niter, ATA, b, x, small_art_kspace1, acq)
r = b - ATA(x);
p = r;
rsold = r' * r;

for idx = 1:length(b)
    Ap = ATA(p);
    alpha = rsold / (p' * Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    rsnew = r' * r;
    if sqrt(rsnew) < 1e-25 || idx>niter
        break;
    end
    p = r + (rsnew / rsold) * p;
    err(idx) = norm(rsold-rsnew);
    rsold = rsnew;
    
end

return;


function [new_kspace] = sense_op_t(coils_ind,x,loc_mask,m,n,no_c)

% put the appr. locations of the conc. kspace with the incoming data
kk1 = zeros(m,n,no_c,'single');
kk1(loc_mask) = x;


% go back to image domain
kk2 = fftshift(ifft2(ifftshift(kk1))); %% JAM change for consitency


% conjugate sens multiplication
for abc = 1:no_c
    ev1(:,:,abc) = conj(coils_ind(:,:,abc)).*kk2(:,:,abc);
end

ev2 = sum(ev1,3);

new_kspace = ev2(:);
return;



function [new_kspace] = sense_op(sensMap, ksp, nRO, nPE, nCh, maskAcquired)
kk1 = reshape(ksp, nRO, nPE);

% sens multiplication
ev1 = sensMap .* repmat(kk1,[1 1 nCh]);

% go back to k-space
ev2 = fftshift(fft2(ifftshift(ev1))); % JAM change for consistency

% taking the acq points from k-space
new_kspace = ev2(maskAcquired);

return;


