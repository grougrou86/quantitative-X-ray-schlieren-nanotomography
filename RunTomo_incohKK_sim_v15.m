%% input
% jsonPath =[];
topPath = 'G:\20240911_BMOL\MATLAB'; % G:\20240521_BMOL\MATLAB
addpath(  genpath(topPath))

%[tH, jsonPathOut] = tomoHandle([],'KK');
%disp(jsonPathOut)

%MAIN_PATH='G:\20240911_BMOL\mixture';name='\NO_FILTER_8374.000';percent_cropp=0.3;
%MAIN_PATH='G:\20240911_BMOL\silicon';name='\NO_FILTER_8374.000';percent_cropp=0.3;
%MAIN_PATH='G:\20240911_BMOL\NCM2';name='\NCM_8374.000';percent_cropp=0.4;
%MAIN_PATH='G:\20240911_BMOL\discharge_sample_2_with_filter';name='\NO_FILTER_8374.000';percent_cropp=0.6;
%MAIN_PATH='G:\20240911_BMOL\charged';name='\NO_FILTER_8374.000';percent_cropp=0.8;
MAIN_PATH='F:\simulation_xray_kk\sim1_noise_10000_phase_1e-05_abs_1e-06';name='\SIM_8333.000';percent_cropp=1;

jsonsetup=['G:\20240911_BMOL\setup_info_with_computed_PSF.json'];

[tH, jsonPathOut] = tomoHandle(jsonsetup,'KK');
disp(jsonPathOut)

%0.7;

%% load cutoff
%titleStr = sprintf('Choose cutoff data file for %d eV',tH.setup.source.energy_eV);
%[f1, p1] = uigetfile('result.mat',titleStr,'../../../');
%cutoff_file = fullfile(p1, f1);

cutoff_file = [MAIN_PATH name '_GSCAN_0'];
fid = fopen([MAIN_PATH '\scan_info.json']);

raw = fread(fid, inf);
str = char(raw');
fclose(fid);
val = jsondecode(str);
%error('crop FOV during cutoff analysis ')
[fourier_position,fourier_intensity,dark_position]=cut_off_analysis_general_v2(cutoff_file,val.start,val.step);
% {
selection=25:70;
fourier_position=fourier_position(selection);
fourier_intensity=fourier_intensity(selection);
%}
p=struct;
p.fourier_position  = fourier_position * 1e3; % um;
p.fourier_intensity = fourier_intensity;  % minus is required as we scanned max to min (diff is negative)
p.wavelength = Etowl(tH.setup.source.energy_eV * 1e-3)* 1e-3; % um

p.cutoff = val.result_cutoff * 1e3; % um


figure,
plot(p.fourier_position, p.fourier_intensity);
xlabel('fourier position (um)')
xline(p.cutoff,'-','cutoff position')

% if you DO NOT want to use gaussian fit for the pattern
%p = rmfield(p,'gaussFit'); warning('gaussian filtering removed');


%% setData
[~,~,ext] = tH.setData('TOMO_00000.tif',[MAIN_PATH name '_00001']);

%%% flyscan angle adjustment
field_location = tH.samFiles(1).folder;

%%% set rotAngleMax
fname = fullfile(field_location, 'metadata.txt');
fid = fopen(fname);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
val = jsondecode(str);

%%% 0-180 deg pairing
rotAnglelib = linspace(val.stop_angle-360, val.stop_angle, length(tH.samFiles));
[minAngleDiff, ind180] = min(abs(rotAnglelib - 180));
ind360 = 2*(ind180 - 1);

% update
tH.samFiles = tH.samFiles(1:ind360);
tH.rotAnglelib = rotAnglelib(1:ind360);
tH.rotAngleMax = rotAnglelib(ind360);

%%% set tH.darkFiles (put dark field image at the jsonPath)
%tH.darkFiles = dir(fullfile(fileparts(jsonPathOut),['*',ext]));
tH.darkFiles = dir(dark_position);
assert(~isempty(tH.darkFiles), 'put dark field image at the ''jsonPath''')

%%% setFOV

FOV = tH.setFOV(0.75);

%%% get a raw Sinogram
%--------------------------------------------------------------------------------------------------------------------------
indSpacing = 1;
%error('change spacing back to 1');
if indSpacing ~= 1
    warning('not all angles used');
end

samInds = 1:indSpacing:(ind180-1);
%error('cahnge back')
samInds = [samInds, samInds + ind180 - 1];
[rawSino,raw_BG] = tH.getSinogram(samInds);
% rawSino = tH.getSinogram();
angleLib = tH.rotAnglelib(samInds);
% angleLib = tH.rotAnglelib();
% figure, contrastVis(rawSino);

%% deconvolution
rawSino = tH.deconvSinogram(rawSino); % BWcrop inside --> tH.setup.imagePixelSize changed
rawSino = gather(rawSino);

raw_BG = tH.deconvSinogram(raw_BG); % BWcrop inside --> tH.setup.imagePixelSize changed
raw_BG = gather(raw_BG);
% figure, contrastVis(rawSino);

%%% bgnorm
[yy, xx, Nang0] = size(rawSino);
bgWingSize = round(tH.FOV.bgAddXFOV * xx / tH.FOV.size(2));
bgMask = false(yy, xx);
bgMask(:, [1:bgWingSize, end-bgWingSize+1:end]) = true;
rawSino = reshape(rawSino, yy*xx, Nang0);
bgVal   = mean( rawSino(bgMask,:) , 1);

rawSino = rawSino./bgVal;
rawSino = reshape(rawSino, yy, xx, Nang0);
% figure, contrastVis(rawSino);


%%% set field recon parameters
p.pixel_size = tH.setup.imagePixelSize * tH.FOV.size(2)/xx * 1e-3; % um
% scintPix =  tH.setup.detector.camera.pixelSize * tH.setup.detector.camera.binning / tH.setup.detector.objectiveLens.magnification;
f_ZP = tH.setup.afterSample.zonePlate.diameter * tH.setup.afterSample.zonePlate.outermostWidth / Etowl(tH.setup.source.energy_eV * 1e-3); %focal length from feature size and diameter
p.ZP_focal_length = f_ZP * 1e-3;  % um
% mag_xray = scintPix /tH.setup.imagePixelSize;
% params_recon.sample_lens_dist = f_ZP * (1+1/mag_xray) * 1e-3; % um
rawSino_0 = rawSino;
%error('stop')
%%
%{
for pp=1:size(rawSino,1)
    display(['strip : ' num2str(pp) ' / ' num2str(size(rawSino,1))])
    rawSino(pp,:,:)=removeStripes('h',gpuArray(squeeze(rawSino_0(pp,:,:))),2,'db45',3);
end
%}
%% axisShift and num iter tuning
offset=-5;%-13;
elipsity_pos=0;
elipsity_strength=0;
possible_val=-25:0.5:25;
rms_val=[];
rms_val2=[];
rms_val3=[];
for offset=possible_val;

    num_itter = 1;
    badFrames = []; % affected by cosmic radiation..

    % + 5* 5;
    %tH.axisOffset = -6 + 0.5 *-2 + 0.2*-2;
    tH.axisOffset = offset;%+ 0.5 *-2 + 0.2*-2;
    % p.cutoff = C{2}(matchInd)*1e3 + 5*0; % um

    %%% bad angle, and its 180 deg counterpart removal
    Tang180 = Nang0/2-1;
    badFrames = mod(badFrames,Tang180);
    badFrames = unique([badFrames, badFrames+Tang180]);
    validAngInd = true(1,Nang0);
    validAngInd( badFrames ) = false;
    if ~isempty(badFrames)
        warning('nonempty badFrames')
    end

    yind = floor(size(rawSino,1)/2);
    if length(possible_val)>1
        rawSino_rec=rawSino(yind,:,validAngInd);
    else
        rawSino_rec=permute(removeStripes('h',gpuArray(squeeze(rawSino(yind,:,validAngInd))),4,'db45',3),[3 1 2 ]);
        num_itter = 15;
    end
    %rawSino_rec=circshift(rawSino,[0 tH.axisOffset 0]);
    rawSino_rec=s_fft(rawSino_rec,2);
    kkx=(1:size(rawSino_rec,2))-floor(size(rawSino_rec,2)/2)-1;kkx=kkx./size(rawSino_rec,2);
    kkx=reshape(kkx,1,[]);
    kky=(1:size(rawSino_rec,3))-floor(size(rawSino_rec,3)/2)-1;kky=kky./size(rawSino_rec,3);
    kky=reshape(kky,1,1,[]);


    %rawSino_rec=rawSino_rec.*exp(-1i.*2.*pi.*kkx.*(offset+elipsity_strength*sin(2*pi*elipsity_pos+4*pi*kky)));

    rawSino_rec=rawSino_rec.*exp(-1i.*2.*pi.*kkx.*(offset));
    if length(possible_val)>1
        rawSino_rec=rawSino_rec.*(1-exp(-40.^2*(kkx.^2)));
    end

    rawSino_rec=(s_ifft(rawSino_rec,2));
    rawSino_rec=...
        rawSino_rec(...
        (1+floor(size(rawSino_rec,1)*((1-percent_cropp)/2))):...
        (size(rawSino_rec,1)-floor(size(rawSino_rec,1)*((1-percent_cropp)/2))),...
        (1+floor(size(rawSino_rec,2)*((1-percent_cropp)/2))):...
        (size(rawSino_rec,2)-floor(size(rawSino_rec,2)*((1-percent_cropp)/2))),:);

    %%% run
    fieldReconInputParams = {num_itter, p, 'figShow', false, 'DPCdata', struct([])};


    [tomoOut, Eout, DPCdata] = sino2tomo(gpuArray(rawSino_rec), angleLib(validAngInd), 0, bgWingSize, fieldReconInputParams, 1,rawSino_rec(:,:,1)*0+1);
    %[tomoOut, ~, DPCdata] = sino2tomo(gpuArray(rawSino(yind,:,validAngInd)), angleLib(validAngInd), tH.axisOffset, bgWingSize, fieldReconInputParams, 1);
    % fieldReconinputParams = {num_itter,p,'figShow',figShow,'DPCdata', DPCdata}
    % figure, imagesc(squeeze(rawSino(yind,:,:)));
    % figure, contrastVis(rawSino(:,:,669));
    rms_val(end+1)=sqrt(mean(imag(tomoOut).^2,'all'));
    rms_val2(end+1)=sqrt(mean(abs(imag(tomoOut)),'all'));
    rms_val3(end+1)=sqrt(mean(real(Eout).^2,'all'));

    rms_val
    %error('pause')
end
%rms_val4=rms_val./rms_val3;
rms_val4=rms_val2;
figure; plot(rms_val4);title('cost function to determine the best shift ')
%drawnow;
%error('stop')
%% tomoRecon
[~,best_id]=min(rms_val4(:));
offset=possible_val(best_id);
%tH.axisOffset =possible_val(best_id);
%display(['Best shift is : ' num2str(tH.axisOffset)])
%figure; plot(possible_val,rms_val4);
rawSino_rec=rawSino;
raw_BG_rec=raw_BG;
%rawSino_rec=circshift(rawSino,[0 tH.axisOffset 0]);
rawSino_rec=s_fft(rawSino_rec,2);
raw_BG_rec=s_fft(raw_BG_rec,2);

kkx=(1:size(rawSino_rec,2))-floor(size(rawSino_rec,2)/2)-1;kkx=kkx./size(rawSino_rec,2);
kkx=reshape(kkx,1,[]);
kky=(1:size(rawSino_rec,3))-floor(size(rawSino_rec,3)/2)-1;kky=kky./size(rawSino_rec,3);
kky=reshape(kky,1,1,[]);

%rawSino_rec=rawSino_rec.*exp(-1i.*2.*pi.*kkx.*(offset+elipsity_strength*sin(2*pi*elipsity_pos+4*pi*kky)));
raw_BG_rec=raw_BG_rec.*exp(-1i.*2.*pi.*kkx.*(offset));
raw_BG_rec=(s_ifft(raw_BG_rec,2));
raw_BG_rec=...
    raw_BG_rec(...
    (1+floor(size(raw_BG_rec,1)*((1-percent_cropp)/2))):...
    (size(raw_BG_rec,1)-floor(size(raw_BG_rec,1)*((1-percent_cropp)/2))),...
    (1+floor(size(raw_BG_rec,2)*((1-percent_cropp)/2))):...
    (size(raw_BG_rec,2)-floor(size(raw_BG_rec,2)*((1-percent_cropp)/2))),:);

rawSino_rec=rawSino_rec.*exp(-1i.*2.*pi.*kkx.*(offset));
rawSino_rec=(s_ifft(rawSino_rec,2));
rawSino_rec=...
    rawSino_rec(...
    (1+floor(size(rawSino_rec,1)*((1-percent_cropp)/2))):...
    (size(rawSino_rec,1)-floor(size(rawSino_rec,1)*((1-percent_cropp)/2))),...
    (1+floor(size(rawSino_rec,2)*((1-percent_cropp)/2))):...
    (size(rawSino_rec,2)-floor(size(rawSino_rec,2)*((1-percent_cropp)/2))),:);

%%
rawSino_rec_out=fine_register_fields_TV_v5(rawSino_rec);

%error('stop');
%%
num_itter = 7;%7;
fieldReconInputParams = {num_itter, p, 'figShow',true,'DPCdata',struct([])};

recon_aligned=true;

crp_final=5;
rawSino_rec_out_final=rawSino_rec_out(1+crp_final:end-crp_final,1+crp_final:end-crp_final,:);
%raw_BG_rec_out_final=raw_BG_rec_out(1+crp_final:end-crp_final,1+crp_final:end-crp_final,:);
if recon_aligned
    [yy, xx, Nang0] = size(rawSino_rec_out_final);
else
    [yy, xx, Nang0] = size(rawSino_rec);
end

opt_fun=@(x) try_RI_recon(x,rawSino_rec_out_final,yy,xx,Nang0,angleLib,bgWingSize,fieldReconInputParams);
a=0;
b=1;%1;
itt_gss=5;%5;
invphi = (sqrt(5) - 1) / 2 ;
comp_c=true;
comp_d=true;
for gss_crr =1:itt_gss
    c = b - (b - a) * invphi;
    d = a + (b - a) * invphi;
    if comp_c
        fc=opt_fun(c);
    end 
    if comp_d
        fd=opt_fun(d);
    end
    if fc < fd
        b = d;
        fd=fc;
        comp_c=true;
        comp_d=false;
    else
        a = c;
        fc=fd;
        comp_c=false;
        comp_d=true;
    end
    display(['a = ' num2str(a) '& b = ' num2str(b)]);
end
optimal_shift=(b + a) / 2;

[optimal_cost,RI_optimal]=opt_fun(optimal_shift);
RI=RI_optimal;

figure; plot(fig_ring_validity_tester(RI_optimal));

figure, contrastVis(real(RI_optimal)); axis image; colorbar; %title(titleStr, 'Interpreter', 'none')
figure, contrastVis(imag(RI_optimal)); axis image; colorbar;

%%
img_window=blackman(size(RI,1)).*blackman(size(RI,2))'.*reshape(blackman(size(RI,3)),1,1,[]);
figure; orthosliceViewer(log(abs(s_fftn(real(RI.*img_window)))+0.0003));
%% filter
%{
sz_grid = ndgrid_matSizeIn(size(filtered_RI),true,'centerZero');

filtered_RI=s_fftn(RI);
percent_vert=0.75;
%filtered_RI=filtered_RI.*(exp(-10.^2.*(exp(-100.^2*(sz_grid{1}).^2)).*((sz_grid{3}).^2+(sz_grid{2}).^2)));
filtered_RI=s_ifftn(filtered_RI);

figure, contrastVis(real(filtered_RI)); axis image; colorbar;
%}
%% save
%{
saveTag = 'v5_herve';
% v4: save rawSino, bgWingSize, validAngInd, and delta = real(RI);
% v3: save fieldReconInputParams;

saveStr = sprintf('tomoResult_%s.mat', saveTag );
saveStr = fullfile( tH.samFiles(1).folder, saveStr );
fprintf('[%s] save ... ', datetime)
save(saveStr,'tH','fieldReconInputParams','RI','rawSino',"bgWingSize","validAngInd","-v7.3","-nocompression");
fprintf('done\n')
%}
%%
%save([MAIN_PATH name '_RI_result.mat'],'RI','rawSino_rec_out_final','-v7.3');

function [err_coo,RI]=try_RI_recon(shift_axis,rawSino_rec_out_final,yy,xx,Nang0,angleLib,bgWingSize,fieldReconInputParams)
refocus_curve=0;%-1;
RI = zeros(yy,xx,xx,'single');
bytesPerSlice = (xx*Nang0) * 8; % 4 bytes for single
gpu=gpuDevice;
maxChunkSize = 500;%floor( gpu.AvailableMemory / bytesPerSlice / 30 ); % last factor is a heuristic safety value.
Npart = ceil( yy / maxChunkSize);
chunks = round( linspace(1, yy+1, Npart+1) );
%shift_axis=0;%0.3;
%error('choose shift automaticaly using the cross correlation method');
for ii = 1:Npart
    fprintf( '[%s] chunk start: %03d/%03d\n', datetime,ii,Npart )
    thisChunk = chunks(ii) : chunks(ii+1)-1;
    nmmm=1;
    validAngInd=[(1:nmmm:floor(size(rawSino_rec_out_final,3)/2)) (floor(size(rawSino_rec_out_final,3)/2)+1):nmmm:size(rawSino_rec_out_final,3)];

    RI(thisChunk,:,:) = sino2tomo(gpuArray(rawSino_rec_out_final(thisChunk,:,validAngInd)), angleLib(validAngInd), shift_axis, bgWingSize, fieldReconInputParams, 0,1,refocus_curve);

end
coo=fig_ring_validity_tester(RI);
coo=coo(1:10);
coo=(coo./mean(coo(:)))-1;
%err_coo=std(coo(:));
err_coo=mean(abs(coo(:)));

%error('automate finding the shift without computing every shift by finding slope and deducing the shift ');
end




