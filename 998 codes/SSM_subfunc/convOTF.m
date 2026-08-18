function [deconv_output, meanSNRMap] = deconvOTF(FrefImg,scintOTF,signalkR,varargin)
% p = inputParser();
% addParameter(p,'iterMax',1000);
% parse(p,varargin{:});
% 
% iterMax = p.Results.iterMax;
meanSNRMap=1;
deconv_output=s_ifft2(FrefImg.*(scintOTF));%.*(abs(scintOTF)>0.2));
%error('modify')

return ;
%%
[yy,xx,zz] = size(FrefImg);
signalWindow = ~mk_ellipse(signalkR*xx,signalkR*yy,xx,yy);
% figure,imagesc(signalWindow),axis image

if all(signalWindow,'all')
    invWienerOTF = 1./scintOTF;
else
    FrefImg = reshape(FrefImg,[],zz);

    Sf = abs(FrefImg(~signalWindow,:)).^2;
    bgMean = mean( Sf, 1);
    bgStd  = std ( Sf,0,1 );
    clear Sf

    FrefImg = reshape(FrefImg,yy,xx,zz);
    bgMean = reshape(bgMean,1,1,zz);
    bgStd  = reshape(bgStd,1,1,zz);

%     SNR = abs(FrefImg).^2./bgMean;
    SNR = (abs(FrefImg).^2 - bgMean)./bgStd/3; % factor of 3 is heuristrical value.
    SNR(SNR<0) = 0;
    %     figure,tomoHandle.show( (1 + 1 ./ (abs(scintOTF).^2.*SNR) ).^-1)
    invWienerOTF = 1./scintOTF .* (1 + 1 ./ (abs(scintOTF).^2.*SNR) ).^-1;
    meanSNRMap = mean(SNR,3);
    clear SNR
    
    invWienerOTF( isnan(invWienerOTF) ) = 0;    
end

invWienerOTF = scintOTF./(scintOTF.^2+0.1);

FrefImg = FrefImg .* invWienerOTF;    % Wiener deconvolution
clear invWienerOTF

% deconv_output = fftshift(ifft2(ifftshift(Ftemp.*signalWindow)));
% imagesc( (1 + 1 ./ (abs(scintOTF).^2.*SNR) ).^-1)
% imagesc( abs(invWienerOTF(:,:,1)) )
% imagesc( abs(SNR(:,:,1)) )

deconv_output = fftshift(ifft2(ifftshift(FrefImg)));
% input = fftshift(ifft2(ifftshift(FrefImg1)));
% 
% subplot(121),imagesc(log10(abs(deconv_output(:,:,1)))),axis image
% subplot(122),imagesc(abs(input)),axis image


% for iter = 1:1:iterMax
%     Fdeconv_output = Ftemp.*signalWindow;
%     
%     % imagesc(log10(abs(Feximg_deconv))),axis image;
%     % imagesc(abs(FI_PSF.*signalWindow)),axis image;
%     deconv_output = fftshift(ifft2(ifftshift(Fdeconv_output)));
%     %     disp(min(eximg_deconv(:)))
%     minval = min(deconv_output(:));
%     deconv_output(deconv_output<0) = 0;
%     if abs(minval) < mean(deconv_output(:))/2^16
%         deconv_output = real(deconv_output);
%         break;
%     end
%     Ftemp = fftshift(fft2(ifftshift(deconv_output)));
% end

