function [p,H] = iradon_filterPart_v2(p, filter, frequencyScaling, spatialDimThatNormalToRotationAxis,SNRplot)
%IRADON Inverse Radon transform.
%   I = iradon(R,THETA) reconstructs the image I from projection data in
%   the 2-D array R.  The columns of R are parallel beam projection data.
%   IRADON assumes that the center of rotation is the center point of the
%   projections, which is defined as ceil(size(R,1)/2).
%
%   THETA describes the angles (in degrees) at which the projections were
%   taken.  It can be either a vector containing the angles or a scalar
%   specifying D_theta, the incremental angle between projections. If THETA
%   is a vector, it must contain angles with equal spacing between them.
%   If THETA is a scalar specifying D_theta, the projections are taken at
%   angles THETA = m * D_theta; m = 0,1,2,...,size(R,2)-1.  If the input is
%   the empty matrix ([]), D_theta defaults to 180/size(R,2).
%
%   IRADON uses the filtered backprojection algorithm to perform the
%   inverse Radon transform.  The filter is designed directly in the
%   frequency domain and then multiplied by the FFT of the projections.
%   The projections are zero-padded to a power of 2 before filtering to
%   prevent spatial domain aliasing and to speed up the FFT.
%
%   I = IRADON(R,THETA,INTERPOLATION,FILTER,FREQUENCY_SCALING,OUTPUT_SIZE)
%   specifies parameters to use in the inverse Radon transform.  You can
%   specify any combination of the last four arguments.  IRADON uses
%   default values for any of these arguments that you omit.
%
%   INTERPOLATION specifies the type of interpolation to use in the
%   backprojection. The default is linear interpolation. Available methods
%   are:
%
%      'nearest' - nearest neighbor interpolation
%      'linear'  - linear interpolation (default)
%      'spline'  - spline interpolation
%      'pchip'   - shape-preserving piecewise cubic interpolation
%
%   FILTER specifies the filter to use for frequency domain filtering.
%   FILTER is a string or a character vector that specifies any of the
%   following standard filters:
%
%   'Ram-Lak'     The cropped Ram-Lak or ramp filter (default).  The
%                 frequency response of this filter is |f|.  Because this
%                 filter is sensitive to noise in the projections, one of
%                 the filters listed below may be preferable.
%   'Shepp-Logan' The Shepp-Logan filter multiplies the Ram-Lak filter by
%                 a sinc function.
%   'Cosine'      The cosine filter multiplies the Ram-Lak filter by a
%                 cosine function.
%   'Hamming'     The Hamming filter multiplies the Ram-Lak filter by a
%                 Hamming window.
%   'Hann'        The Hann filter multiplies the Ram-Lak filter by a
%                 Hann window.
%   'none'        No filtering is performed.
%
%   FREQUENCY_SCALING is a scalar in the range (0,1] that modifies the
%   filter by rescaling its frequency axis.  The default is 1.  If
%   FREQUENCY_SCALING is less than 1, the filter is compressed to fit into
%   the frequency range [0,FREQUENCY_SCALING], in normalized frequencies;
%   all frequencies above FREQUENCY_SCALING are set to 0.
%
%   OUTPUT_SIZE is a scalar that specifies the number of rows and columns
%   in the reconstructed image.  If OUTPUT_SIZE is not specified, the size
%   is determined from the length of the projections:
%
%       OUTPUT_SIZE = 2*floor(size(R,1)/(2*sqrt(2)))
%
%   If you specify OUTPUT_SIZE, IRADON reconstructs a smaller or larger
%   portion of the image, but does not change the scaling of the data.
%
%   If the projections were calculated with the RADON function, the
%   reconstructed image may not be the same size as the original image.
%
%   [I,H] = iradon(...) returns the frequency response of the filter in the
%   vector H.
%
%   Class Support
%   -------------
%   R can be double or single. All other numeric input arguments must be
%   double. I has the same class as R. H is double.
%
%   References
%   -----------
%   [1] A. C. Kak, Malcolm Slaney, "Principles of Computerized Tomographic
%   Imaging", IEEE Press 1988.
%
%   Example 1:
%   ----------
%   % Compare filtered and unfiltered backprojection.
%
%       P = phantom(128);
%       R = radon(P,0:179);
%       I1 = iradon(R,0:179);
%       I2 = iradon(R,0:179,'linear','none');
%       subplot(1,3,1), imshow(P), title('Original')
%       subplot(1,3,2), imshow(I1), title('Filtered backprojection')
%       subplot(1,3,3), imshow(I2,[]), title('Unfiltered backprojection')
%
%   Example 2:
%   ----------
%   % Compute the backprojection of a single projection vector. The IRADON
%   % syntax does not allow you to do this directly, because if THETA is a
%   % scalar it is treated as an increment.  You can accomplish the task by
%   % passing in two copies of the projection vector and then dividing the
%   % result by 2.
%
%       P = phantom(128);
%       R = radon(P,0:179);
%       r45 = R(:,46);
%       I = iradon([r45 r45], [45 45])/2;
%       imshow(I, [])
%       title('Backprojection from the 45-degree projection')
%
%   See also FAN2PARA, FANBEAM, IFANBEAM, PARA2FAN, PHANTOM, RADON.

%   Copyright 1993-2020 The MathWorks, Inc.

%%%
if nargin< 5
    SNRplot = [];
end

%%% make the fft dimension front
realTF = isreal(p);
permuteVec = 1:ndims(p);
permuteVec = circshift(permuteVec, -spatialDimThatNormalToRotationAxis + 1);
p = permute( p, permuteVec);

[yy,xx,zz] = size(p);
H = designFilter(filter, yy, frequencyScaling);
% H = designFilter('hamming', yy, xWindowSize);
% plot(H)

if strcmpi(filter, 'none')
    return;
end

%%% padding and fft
yy2 = length(H);
% p = padarray(p,yy2 - yy , 0 , 'post'); %% padding
p = cat(1, p, flip(p,1) );
% p = padarray(p,yy2 - yy , 0 , 'post'); %% padding

p = fft(p,[],1);               % p holds fft of projections

% imagesc(log10(abs( p(:,:,1) ))), axis image
% imagesc(abs( p(:,:,1) )), axis image


%%% Wiener deconvolution
if strcmpi(filter, 'wiener')
    bgVec   = false(yy2,1);
    bgRange = [frequencyScaling/2, 1-frequencyScaling/2]*(yy2-1) + 1; 
    bgRange = [floor(bgRange(1)), ceil(bgRange(2))];
    bgVec(bgRange(1):bgRange(2)) = true;
    
    if isempty(SNRplot)
        Sf = abs(p).^2;
        bgMean = mean( Sf(bgVec,:,:),  'all');
        bgStd  = std(  Sf(bgVec,:,:),0,'all');
        % imagesc(Sf(:,:,100)-bgVal)
        Sf = (Sf-bgMean)/bgStd; % Sf->= SNR
        %     Sf = (Sf-bgMean)/bgMean; % Sf->= SNR
        Sf(Sf<0) = 0;
        % figure,imagesc(log10(abs(Sf)))
        H = H.*(1 + abs(H).^2 ./ Sf ).^-1;
        % figure,imagesc(H)
    else
        %%% interp SNR_kx
        N = length(SNRplot);
        f0 = ((1:N)-mcoor(N))/N;

        N = length(H);
        f1 = [0:N-mcoor(N),1-mcoor(N):-1]/N;

        Sf = interp1(f0,SNRplot,f1,'linear');
        Sf( Sf<0 | isnan(Sf) ) = 0;
        % plot(f0,log10(SNRplot),f1,log10(Sf))
        % plot(log10(Sf))

        Sf = reshape(Sf,size(H,1),size(H,2));
        Sf = Sf/10;
        H = H.*(1 + abs(H).^2 ./ Sf ).^-1;
        % plot(H)
    end
    

% imagesc( H(:,:,160) )
% imagesc( (1 + abs(H).^2 ./ Sf ).^-1 )
% imagesc( H.*(1 + abs(H).^2 ./ Sf(:,:,160) ).^-1 )
% plot(H)
end

%     bgRange = (1/2 + [-1/yy, 1/yy])*yy2;
%     bgRange = [floor(bgRange(1)), ceil(bgRange(2))];
%     Sf = mean( abs(p).^2, [2,3] );
%     bgMean = mean( Sf( bgRange(1):bgRange(2)) );
%     bgStd  = std(  abs(p(bgRange(1):bgRange(2),:,:)).^2,0,'all');
%     % imagesc(Sf(:,:,100)-bgVal)
%     Sf = (Sf-bgMean)/bgStd; % Sf->= SNR
% %     Sf = (Sf-bgMean)/bgMean; % Sf->= SNR
%     Sf(Sf<0) = 0;
%     % imagesc(Sf(:,:,160))
%     H = H.*(1 + abs(H).^2 ./ Sf ).^-1;
    

%%% deconv
H(isnan(H)) = 0;
p = p.* H;
% imagesc( log10(abs(p)))

% p = bsxfun(@times, p, H); % faster than for-loop
if realTF
    p = ifft(p,[],1,'symmetric');  % p is the filtered projections
else
    p = ifft(p);  % p is the filtered projections
end

p = reshape(p,length(H) ,[]);
p(yy+1:end,:) = [];      % Truncate the filtered projections
p = reshape(p,yy,xx,zz);

p = ipermute( p, permuteVec);
% imagesc( squeeze(p))

%----------------------------------------------------------------------

%======================================================================
function filt = designFilter(filter, len, d)
% Returns the Fourier Transform of the filter which will be
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections


% order = max(64,2^nextpow2(2*len));
order = 2*len;

if strcmpi(filter, 'none')
    filt = ones(1, order);
    return;
end

% First create a bandlimited ramp filter (Eqn. 61 Chapter 3, Kak and
% Slaney) - go up to the next highest power of 2.

n = 0:(order/2); % 'order' is always even. 
filtImpResp = zeros(1,(order/2)+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/4; % Set the DC term 
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2); % Set the values for odd n
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 

filt = 2*real(fft(filtImpResp)); 
filt = filt(1:(order/2)+1);

w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

switch filter
    case {'ram-lak','wiener'}
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter
