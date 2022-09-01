% PJB testing an SMS demo

% Load the brain 
load mri
sliceA = double(squeeze(D(:,:,1,6)));
sliceB = double(squeeze(D(:,:,1,21)));

% Normalize
sliceA = sliceA ./ max(double(D(:)));
sliceB = sliceB ./ max(double(D(:)));

figure(1)
subplot(1,2,1)
imagesc(sliceA, [0, 1])
title('True slice A');
subplot(1,2,2)
imagesc(sliceB, [0, 1])
title('True slice B');

% Assume 2 coils, with sensitivity profiles C1 and C2
% Synthesize some sens profiles
N = 128;
if true
    % Spatially varying sensitivity profile
    [YY, XX] = meshgrid(1:N, 1:N);
    C1A = sin(pi/2*YY./N) .* sin(pi/2*XX./N);
    % Rotate coil for all the others
    C1B = 0.8 * C1A(:, end:-1:1);
    C2A = 0.6 * C1A(end:-1:1, :);
    C2B = 1.2 * C1A(end:-1:1, end:-1:1);
else
    % Uniform sensitivity profile
    C1A = ones(N, N);
    % Scale coil for all the others
    C1B = 0.8 * C1A;
    C2A = 0.6 * C1A;
    C2B = 1.2 * C1A;
end



% Display coil sensitivity profiles
figure(2)
subplot(2,2,1)
imagesc(C1A, [0,1])
title('Coil 1, measuring slice A');
subplot(2,2,2)
imagesc(C1B, [0,1])
title('Coil 1, measuring slice B');
subplot(2,2,3)
imagesc(C2A, [0,1])
title('Coil 2, measuring slice A');
subplot(2,2,4)
imagesc(C2B, [0,1])
title('Coil 1, measuring slice B');

% In image domain, the images we measure are the product of the both slices
% and the sensitivity profiles
S1 = C1A.*sliceA + C1B.*sliceB;
S2 = C2A.*sliceA + C2B.*sliceB;

if true % add noise?
    noiseval = .02;
    S1 = S1 + randn(N, N).*noiseval;
    S2 = S2 + randn(N, N).*noiseval;
end

% Show measured data
figure(3)
subplot(1,2,1)
imagesc(S1)
title('Measured coil 1');
subplot(1,2,2)
imagesc(S2)
title('Measured coil 2');


% Point by point, we setup a matrix and solve
% Note - much more efficient to set up and solve it all at once
IA = zeros(N, N);
IB = zeros(N, N);
for idx=1:N
    fprintf('%d/%d\n', idx, N);
    for jdx=1:N
        
        % Solve a matrix of 2 eqns 2 unknowns, format Ax = b
        % A is sesnsitivity encoding C
        % x = is the true signal 
        % b = s the measured signal 
        % Use matlab's Ax=b solver: x = A\b, or image = C\S
        
        clear S C;
        S(1) = S1(idx, jdx);
        S(2) = S2(idx, jdx);
        %if max(S(:)) ~= 0 % Skip the case of both zeros
            
            C(1,1) = C1A(idx, jdx);
            C(1,2) = C2A(idx, jdx);
            C(2,1) = C1B(idx, jdx);
            C(2,2) = C2B(idx, jdx);
            
            tmp = C.'\S.';
            IA(idx, jdx) = tmp(1);
            IB(idx, jdx) = tmp(2);
        %end
        
    end
end

% Show reconstructed iamges
figure(4)
subplot(1,2,1)
imagesc(IA, [0, 1])
title('Reconstructed slice A');
subplot(1,2,2)
imagesc(IB, [0, 1])
title('Reconstructed slice B');





