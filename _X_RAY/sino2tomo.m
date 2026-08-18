function [tomoOut,Eout, DPCdata] = sino2tomo(rawSino, angleLib, axisShift, bgWingSize, fieldReconinputParams, figShow,BG,refocus_curve,remove_ring)
% fieldReconinputParams = {num_itter,p,'figShow',figShow,'DPCdata', DPCdata}
if nargin < 6
    figShow = false;
end

Nang = length(angleLib);

if size(rawSino,3) == Nang
    thetaInd1 = 1:Nang/2;
    thetaInd2 = Nang/2+1:Nang;

    sino1 = circshift_KR(rawSino, [0,axisShift,0]);
    sino2 = flip(sino1(:,:,thetaInd2),2);
    sino1 = sino1(:,:,thetaInd1);
    clear rawSino
    % contrastVis(sino1 + sino2);
    % Eout0 = field_recon_routine_v15(sino1, sino2, 0, num_itter, params_recon, false);
    shift=0;%0.3;
    %error('use the TV algorithm and autodiff to get the min fast');
    %sino1=sino1(:,:,1);sino2=sino2(:,:,1);
    [Esino, DPCdata] = field_recon_routine_v21(real(fourier_shift(sino1,0,shift/2,0)),real(fourier_shift(sino2,0,-shift/2,0)),BG, fieldReconinputParams{:},'figShow',figShow);
    
    %error('stop')
    %[Esino, DPCdata] = field_recon_routine_v18(sino1,sino2,BG, fieldReconinputParams{:});
    
elseif size(rawSino,3) == Nang/2 % reconstructed field input
    Esino = rawSino;
    thetaInd1 = 1:Nang/2;
    thetaInd2 = Nang/2+1:Nang;
    DPCdata = fieldReconinputParams{6};
else
    error('invalid input')
end
% figure,contrastVis(real(Eout0)); 
% contrastVis(imag(Eout0)); 

%%
%{
%filter for resolution  assessement
Esino=cat(3,Esino,flip(Esino,3));
Esino=(s_ifft((1-exp( ...
    -reshape(((1:size(Esino,3))-floor(size(Esino,3)/2)-1)./size(Esino,3),1,1,[]) ...
    .^2*100^2)).*s_fft(Esino,3),3));
Esino=Esino(:,:,1:end/2);
%}
if nargout >= 2
    Eout = Esino; % 
end

%%% bgnorm
% bgVal = mean(Eout0(:, [1:bgWingSize, end-bgWingSize+1:end], :), 2);
% Eout0 = Eout0 - bgVal;

Esino = gpuArray(Esino);
[yy,xx] = size(Esino, [1,2]);

%%% def. bgmask
bgMask = zeros(1,xx);
bgMask([1:bgWingSize, end-bgWingSize+1:end]) = 1;
bgMask = circshift_KR(bgMask, [1, axisShift ,1]);
bgMask = abs(bgMask) > 0.5;

%%% norm
bgVal  = mean(Esino(:, bgMask, : ), 2);
Esino = Esino - bgVal;

%% refocus

% {
    coo=make_coo(size(Esino));
    mult_lin=1e-3;
    aberration = exp(1i.*2.*pi.*(coo{1}.^2+coo{2}.^2).*refocus_curve);
    %aberration;
    Esino=log(s_ifft2(s_fft2(exp(Esino.*mult_lin)).*aberration))./mult_lin;% the use of mult_lin here is to avoid having to do unwrapping while still correcting for refocus artifact but only in a linear way to correct border artifactics
    %Esino=log(s_ifft2(s_fft2(exp(Esino)).*aberration));% original function

    %}

%%% filter
%{
Esino = iradon_filterPart_v2(Esino, "ram-lak", 1, 2);
bgVal  = mean(Esino(:, bgMask, : ), 2);
Esino = Esino - bgVal;
%}
%%% tomo
% yselect = 203;
% axisShift2 = +2;
% test = circshift_KR( bField0(yselect,:,:), [0,axisShift2,0] );

tomoOut  = zeros(yy,xx,xx,'like',Esino(1));
ang180 = angleLib(thetaInd1);

% for yi = 1:yy
%     tomoOut(yi,:,:) = ...
%         -1  * iradon( squeeze(real(Esino(yi,:,:))), ang180, "linear", 'Ram-Lak', 1, xx )...
%         -1i * iradon( squeeze(imag(Esino(yi,:,:))), ang180, "linear", 'Ram-Lak', 1, xx);
%     fprintf("iradon: %05d/%05d\n",yi,yy)
% end
if remove_ring
    decNum=3;5;
    sigma=4.8;2.4;
    wave_name="sym5";
    Esino=...
        1i.*xRemoveStripesVertical(imag(Esino),decNum,wave_name,sigma)...
        +1.*xRemoveStripesVertical(real(Esino),decNum,wave_name,sigma);
end

fprintf("[%s] iradon %05d slices ... ",datetime,yy)
for yi = 1:yy
    %Esino(yi,:,:)=permute(removeStripes('h',squeeze(Esino(yi,:,:)),4,'db25',1),[3 1 2]);
    
    %tomoOut(yi,:,:) = ...
    %    -1i  * iradon( squeeze(real(Esino(yi,:,:))), ang180, "linear", 'None', 1, xx )...
    %    -1 * iradon( squeeze(imag(Esino(yi,:,:))), ang180, "linear", 'None', 1, xx);
    
    
    tomoOut(yi,:,:) = ...
        -1i  * iradon( squeeze(real(Esino(yi,:,:))), ang180,xx)...
        -1 * iradon( squeeze(imag(Esino(yi,:,:))), ang180,xx);
end
fprintf("done\n")

%radius=(((1:xx)-floor(xx/2)-1).^2+((1:xx)-floor(xx/2)-1)'.^2);radius=radius./

%%% bgnorm
samWindowR = floor( (xx-1)/2 - ceil(abs(axisShift)) );
tomoMask = ~mk_ellipse(samWindowR,samWindowR,xx,xx); 

%%% make bg = 0
samWindoWithoutWingR = samWindowR - bgWingSize;
bgMask  = mk_ellipse(samWindoWithoutWingR,samWindoWithoutWingR,xx,xx) & tomoMask;
%             figure,imagesc(bgMask),axis image
tomoOut = reshape(tomoOut, yy, []);
bgVal   = mean( tomoOut(:,bgMask) , 2);
tomoOut = tomoOut - bgVal; % make bg zero for every y-slice
tomoOut = reshape(tomoOut, yy, xx, xx);


%%% physical value
p = fieldReconinputParams{2};
tomoOut  = tomoOut  * (p.wavelength/2/pi /p.pixel_size);

if figShow
    if yy == 1
        titleStr = sprintf("axisShift=%+.2f, num_itter=%d", axisShift, fieldReconinputParams{1});
        subplot(211), contrastVis(real(tomoOut)); axis image; colorbar; title(titleStr, 'Interpreter', 'none')
        subplot(212), contrastVis(imag(tomoOut)); axis image; colorbar;
    else
        figure, contrastVis(real(tomoOut)); axis image; colorbar; %title(titleStr, 'Interpreter', 'none')
        figure, contrastVis(imag(tomoOut)); axis image; colorbar;
    end

    colormap gray;
    set(gcf,'Color', 'w');
    drawnow;
end
