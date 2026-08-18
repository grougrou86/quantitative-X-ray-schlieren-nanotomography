function [goodimg,centerkvec,gof] = PhiShiftKR(Eimg,varargin)
% unit of centerkvec = Fourier space pixelsize.

p = inputParser();
addParameter(p,'bgMask',[]);   
addParameter(p,'useBWlabel', false); 
addParameter(p,'fitWeight',[]); 
addParameter(p,'centerK',[],@(v) isvector(v) || isempty(v));     % if we know the preset center k vector

parse(p,varargin{:});
bgMask    = p.Results.bgMask;
useBWlabel    = p.Results.useBWlabel;
centerK   = p.Results.centerK;
fitWeight = p.Results.fitWeight;

%% dimension check
if isvector(Eimg)
    mfft  = @(x) fftshift(fft(ifftshift(x)));
    mifft = @(x) fftshift(ifft(ifftshift(x)));
    
    if ~iscolumn(Eimg)
        Eimg = Eimg.';
        rowVectorInputTF = true;
    else
        rowVectorInputTF = false;
    end
    
    if ~iscolumn(fitWeight)
        fitWeight = fitWeight.';
    end
    
elseif ismatrix(Eimg)
    mfft  = @(x) fftshift(fft2(ifftshift(x)));
    mifft = @(x) fftshift(ifft2(ifftshift(x)));
    rowVectorInputTF = false;
else
    error('Input must be a vector or a matrix.')
end
[yy,xx]=size(Eimg);

%% BGmask
if isempty(bgMask)
    bgMask = true(yy,xx);
end

if  useBWlabel
    %%% bgMask connectivity test
    [labeledBGMask, NLabel] = bwlabel(bgMask,4);
    if NLabel > 1 % if the mask is composed of islands
        pixelNums = zeros(NLabel,1);
        for ll = 1:NLabel
            pixelNums(ll) = sum(labeledBGMask == ll,'all');
        end
        %     [~,maxInd] = max(pixelNums);
        %
        %     %%% set bgMaskForFFT for the largest island
        %     bgMaskForFFT = (labeledBGMask == maxInd);
        % else
        %     bgMaskForFFT = bgMask;
    else
        pixelNums= 1;
    end
else
    NLabel    = 1;
    pixelNums = 1;
    labeledBGMask = bgMask;    
end

%% fft-based ramp correction (pixel-based)
shiftVec = zeros(NLabel,2);
for nn = 1:NLabel
    bgMaskForFFT = (labeledBGMask == nn);
    % Eimg = Eimg .* bgMask;
    Fimg = mfft(Eimg.* bgMaskForFFT);

    % testMask = bgMask;
    % testMask(:,1:350) = false;
    % Fimg = mfft(Eimg.* testMask);
    % imagesc(angle(Eimg)),axis image
    % imagesc((angle(Eimg.*bgMask))),axis image

    % Fimg = mfft(bgMask);
    % Fimg = mfft(Eimg);
    %   imagesc(abs(Fimg)),axis image

    % imagesc(angle(Eimg.* bgMask))
    % Fmask = ~mk_ellipse(peakRadius,peakRadius,xx,yy);
    if isempty(centerK)
        [~,maxInd] = max(abs(Fimg(:)));
        [dy,dx]=ind2sub([yy,xx],maxInd);
        dky = dy-(floor(yy/2)+1);
        dkx = dx-(floor(xx/2)+1);
    else
        dky = centerK(1)-(floor(yy/2)+1);
        dkx = centerK(2)-(floor(xx/2)+1);
    end

    shiftVec(nn,:) = [-dky,-dkx];
    %                 imagesc(log10(abs(Fimg))),axis image
    %                 imagesc((angle(Eimg_rough))),axis image
    %                 imagesc((angle(Eimg))),axis image
end
meanShiftVec = sum(shiftVec .* pixelNums,1) ./ sum(pixelNums);

% if any(shiftVec(1,:) ~= shiftVec(2,:))
%     disp('different')
% end

[XX,YY] = meshgrid(1:xx,1:yy);
XX = XX - (floor(xx/2)+1);
YY = YY - (floor(yy/2)+1);

% Fimg = circshift( mfft(Eimg), [-dky,-dkx]);
% Eimg_rough = mifft(Fimg);
Eimg_rough = Eimg.*exp(+1i*2*pi*( meanShiftVec(2)*XX/xx + meanShiftVec(1)*YY/yy ));

%% Subpixel registration method
%{
% Foroosh, Hassan, Josiane B. Zerubia, and Marc Berthod. "Extension of
% phase correlation to subpixel registration." IEEE transactions on image processing 11.3 (2002): 188-200.  
% doi: 10.1109/83.988953
% ******* becareful when you use this method with cropped image (which does not fill the entire matrix).

% adjacent pixel selection
direcSelection=[dy-1,dy+1];
[~,choice_pos]=max([Fimg(dy-1,dx),Fimg(dy+1,dx)]);
dy1 =direcSelection(choice_pos);

direcSelection=[dx-1,dx+1];
[~,choice_pos]=max([Fimg(dy,dx-1),Fimg(dy,dx+1)]);
dx1 =direcSelection(choice_pos);

% according to Eq. (22) of Reference.
subDy = [1,1].*Fimg(dy1,dx)./( Fimg(dy1,dx) + [1,-1].*Fimg(dy,dx) );
subDx = [1,1].*Fimg(dy,dx1)./( Fimg(dy,dx1) + [1,-1].*Fimg(dy,dx) );

% imaginary part is noise, and should be (1,-1)
subDy = real(subDy); subDy(abs(subDy)>1)=0;
subDx = real(subDx); subDx(abs(subDx)>1)=0;

%has same sign as [dy1-dy,dx1-dx)
subDy = sign(dy1-dy) * max(subDy);
subDx = sign(dx1-dx) * max(subDx);

dky = dky + subDy;
dkx = dkx + subDx;


[XX,YY] = meshgrid(1:xx,1:yy);
XX = XX - (floor(xx/2)+1);
YY = YY - (floor(yy/2)+1);

% Fimg = circshift( mfft(Eimg), [-dky,-dkx]);
% Eimg_rough = mifft(Fimg);
Eimg_rough = Eimg.*exp(-1i*2*pi*( dkx*XX/xx + dky*YY/yy ));
%                     imagesc(angle(goodimg)),axis image
%}


%% Subpixel registration with poly11 fitting
meanPhase = angle(mean(Eimg_rough(bgMask)));
angtemp = angle(Eimg_rough*exp(-1i*meanPhase));
%                 imagesc(angtemp),axis image
    % imagesc(angtemp.*bgMask),axis image

%%% fitWeight check
if isempty(fitWeight)
elseif all(size(fitWeight) == size(Eimg))
    fitWeight = fitWeight(bgMask);
elseif length(fitWeight) == sum(bgMask)
else
    error('fitWeight error.')
end
fitWeight    = cast(fitWeight,'double');

%%% subpixel fit
% [XX,YY] = meshgrid(1:xx,1:yy);
% XX = XX - (floor(xx/2)+1);
% YY = YY - (floor(yy/2)+1);

if isvector(Eimg)
    [sf,gof] = fit(YY(bgMask),double( angtemp(bgMask) ),'poly1','Weights',fitWeight);
    % plot(sf,YY(bgMask),angtemp(bgMask))
    
    mdytemp = sf.p1;    
    mdxtemp = 0;

elseif ismatrix(Eimg)
    %     [sf,gof] = fit([XX(bgMask),YY(bgMask)],double( gather(angtemp(bgMask)) ),'poly11','Weights',fitWeight);
    %     % plot(sf,[XX(bgMask),YY(bgMask)],angtemp(bgMask))
    %
    %     mdytemp = sf.p01;
    %     mdxtemp = sf.p10;

    %%% Weighted_least_squares
    A = [XX(bgMask), YY(bgMask)]; 
    A = [A,ones(size(A,1),1)];% [N x 3]
    b = angtemp(bgMask); % [N x 1]
    if ~isempty(fitWeight)
        w = sqrt(fitWeight); % [N x 1]
        A = A .* w;
        b = b .* w;
    end

    vecOut = A\b; % solution of Ax = b;
    mdxtemp = vecOut(1);
    mdytemp = vecOut(2);     
end
goodimg = abs(Eimg).*exp(1i.*(angtemp - (mdxtemp*XX+mdytemp*YY)+meanPhase));

% cphi = -mdytemp*YY+mdxtemp*XX;
%                     imagesc(angle(goodimg)),axis image

%%% ouput shift vector based on Fourier-domain pixel 
dky = dky + mdytemp*yy/2/pi;
dkx = dkx + mdxtemp*xx/2/pi;

centerkvec = [dky,dkx];

if rowVectorInputTF
    goodimg    = goodimg.';
    centerkvec = [centerkvec(2),centerkvec(1)];
end






% Ftest = circshift_KR(Fimg,-kvec);
% Etest = fftshift(ifft2(ifftshift(Ftest))); 
% imagesc(angle(Etest))


% BGmasky = logical(bgMask(2:end,:).*circshift(bgMask(1:end-1,:),[0,1]));
% BGmaskx = logical(bgMask(:,2:end).*circshift(bgMask(:,1:end-1),[1,0]));
% 
% %%
% 
% iterthres = 10^-10;
% kmax = 50;
% switch nargin
%     case 1
%         [yy,xx] = size(varargin{1});
%         [goodimg,mdx,mdy] = PhiShiftKR(varargin{1},10,ones(yy,xx));
%         
%     case 2
%         [yy,xx] = size(varargin{1});
%         [goodimg,mdx,mdy] = PhiShiftKR(varargin{1},varargin{2},ones(yy,xx));
%         
%     case 3
%         angimg=varargin{1};
%         peakRadius=varargin{2};
%         bgMask=varargin{3};
%         
% %         if (peakRadius > 1) || (peakRadius < 0)
% %             error('0 <= peakRadius <= 1')
% %         end
%         
%         [yy,xx]=size(angimg);
%         Eimg = exp(1i*angimg);
%         Fimg = fftshift(fft2(ifftshift(Eimg)));
%         Fmask = ~mk_ellipse(peakRadius,peakRadius,xx,yy);
%         [~,mind] = max(abs(Fimg(:)));
%         [dy,dx]=ind2sub([yy,xx],mind);
%         
%         dky = dy-(floor(yy/2)+1);
%         dkx = dx-(floor(xx/2)+1);
%         
%         mdy=2*pi/yy*dky;
%         mdx=2*pi/xx*dkx;
%         
%         Fimg = Fmask .* circshift(Fimg,[-dky,-dkx]);
%         %                 imagesc(log10(abs(Fimg))),axis image
%         Eimg = fftshift(ifft2(ifftshift(Fimg)));
%         angtemp = angle(Eimg);
%         %                 imagesc(angtemp),axis image
%         [XX,YY] = meshgrid(1:xx,1:yy);
%         XX = XX - (floor(xx/2)+1);
%         YY = (floor(yy/2)+1) - YY;
%         
%         BGmasky = logical(bgMask(2:end,:).*circshift(bgMask(1:end-1,:),[0,1]));
%         BGmaskx = logical(bgMask(:,2:end).*circshift(bgMask(:,1:end-1),[1,0]));
%         
%         for kk =1:1:kmax
%             
%             dy_angimg=angle(exp(1i*diff(angtemp,1,1))) .* BGmasky;
%             dx_angimg=angle(exp(1i*diff(angtemp,1,2))) .* BGmaskx;
%             
%             mdytemp = sum(dy_angimg(:))./sum(BGmasky(:));
%             mdxtemp = sum(dx_angimg(:))./sum(BGmaskx(:));
%             
%             if (abs(mdytemp) < iterthres) && (abs(mdxtemp) < iterthres)
%                 break;
%             end
%             
%             cphi = -mdytemp*YY+mdxtemp*XX;
%             angtemp = angle(exp(1i.*(angtemp - cphi)));
%             %                     imagesc(angtemp),axis image
%             mdy = mdy + mdytemp;
%             mdx = mdx + mdxtemp;
%         end
%         
%         if kk == kmax
%             goodimg  = angtemp;
%             mdx = -1;
%             mdy = -1;
%             disp('***** Caution: Iteration does NOT converge: Set proper BGmask *****')
%         else
%             cphi = -mdy*YY+mdx*XX;
%             goodimg = angle(exp(1i.*(angimg - cphi)));
%             mdy = -mdy;
%             %             disp([' PhiShiftKR_iter = ',num2str(kk)])
%         end
%         %                 imagesc(goodimg),axis image
% end