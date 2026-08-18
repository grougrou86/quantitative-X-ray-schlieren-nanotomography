%% input
% jsonPath =[];
topPath = 'G:\20240911_BMOL\MATLAB'; % G:\20240521_BMOL\MATLAB
addpath(  genpath(topPath))

%[tH, jsonPathOut] = tomoHandle([],'KK');
%disp(jsonPathOut)

%MAIN_PATH='G:\20240911_BMOL\mixture';name='\NO_FILTER_8374.000';percent_cropp=0.3;
%MAIN_PATH='G:\20240911_BMOL\silicon';name='\NO_FILTER_8374.000';percent_cropp=0.3;
MAIN_PATH='G:\20240911_BMOL\NCM2';name='\NCM_8374.000';percent_cropp=0.4;
%MAIN_PATH='G:\20240911_BMOL\charged';name='\NO_FILTER_8374.000';percent_cropp=0.8;


files = dir([MAIN_PATH '\*_00*']);
dirFlags = [files.isdir];subFolders = files(dirFlags);subFolderNames = {subFolders(3:end).name} ;
subFolderNames=subFolderNames(6);
%%
for nm=subFolderNames;
    name=['\' nm{1}];
    name=name(1:end-6);

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
    indSpacing = 1;

    samInds = 1:indSpacing:(ind180-1);
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

    %% axisShift and num iter tuning
    offset=-5;%-13;
    elipsity_pos=0;
    elipsity_strength=0;
    possible_val=-50:50;
    rms_val=[];
    rms_val2=[];
    rms_val3=[];
    for offset=possible_val;

    num_itter = 10;
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

    rawSino_rec=rawSino(yind,:,validAngInd);

    %rawSino_rec=circshift(rawSino,[0 tH.axisOffset 0]);
    rawSino_rec=s_fft(rawSino_rec,2);
    kkx=(1:size(rawSino_rec,2))-floor(size(rawSino_rec,2)/2)-1;kkx=kkx./size(rawSino_rec,2);
    kkx=reshape(kkx,1,[]);
    kky=(1:size(rawSino_rec,3))-floor(size(rawSino_rec,3)/2)-1;kky=kky./size(rawSino_rec,3);
    kky=reshape(kky,1,1,[]);

    %rawSino_rec=rawSino_rec.*exp(-1i.*2.*pi.*kkx.*(offset+elipsity_strength*sin(2*pi*elipsity_pos+4*pi*kky)));
    rawSino_rec=rawSino_rec.*exp(-1i.*2.*pi.*kkx.*(offset));
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
    rms_val4=rms_val./rms_val3;
    figure; plot(rms_val4);title('cost function to determine the best shift ')
    %drawnow;
    error('stop')
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
    scann1=-3:1:3;
    scann2=-0.5:0.33/2:0.5;
    rot1=-2:0.5:2;
    rot2=-0.5:0.1:0.5;


    crp_for_align=5;%5;%3.5;
    align_area=rawSino_rec(1+floor(end/crp_for_align):end-floor(end/crp_for_align),1+floor(end/crp_for_align):end-floor(end/crp_for_align),:);
    rawSino_rec_out=rawSino_rec;
    raw_BG_rec_out=raw_BG_rec;

    thetaInd1 = 1:size(align_area,3)/2;
    thetaInd2 = size(align_area,3)/2+1:size(align_area,3);

    align_area1 = align_area;
    align_area2 = gpuArray(single(flip(align_area1(:,:,thetaInd2),2)));
    align_area1 = gpuArray(single(align_area1(:,:,thetaInd1)));

    align_area2_out=(align_area2) ;
    align_area1_out=(align_area1) ;

    coo1=(1:size(align_area1,1))-floor(size(align_area1,1)/2)-1;coo1=coo1./size(align_area1,1);
    coo1=gpuArray(single(reshape(coo1,[],1)));
    coo2=(1:size(align_area1,2))-floor(size(align_area1,2)/2)-1;coo2=coo2./size(align_area1,2);
    coo2=gpuArray(single(reshape(coo2,1,[])));

    coo1_big=(1:size(rawSino_rec,1))-floor(size(rawSino_rec,1)/2)-1;coo1_big=coo1_big./size(rawSino_rec,1);
    coo1_big=gpuArray(single(reshape(coo1_big,[],1)));
    coo2_big=(1:size(rawSino_rec,2))-floor(size(rawSino_rec,2)/2)-1;coo2_big=coo2_big./size(rawSino_rec,2);
    coo2_big=gpuArray(single(reshape(coo2_big,1,[])));

    align_area2 =real(s_ifft2(s_fft2(align_area2).*(1-exp(-(30*coo2).^2))));
    align_area1 =real(s_ifft2(s_fft2(align_area1).*(1-exp(-(30*coo2).^2))));

    filter_border=blackman(size(align_area1,1)).*blackman(size(align_area1,2))';

    ii_opt=0;
    jj_opt=0;
    ann_opt=0;
    figure;
    for aa=1:size(align_area1,3)
        before=align_area1_out(:,:,aa)+align_area2_out(:,:,aa);
        ii_last=0;
        jj_last=0;
        ann_last=0;
        for pp=1:2
            if pp==1
                scann=scann1;
                rot=rot1;
            else
                scann=scann2;
                rot=rot2;
            end
            cost_min=inf;
            for ann=rot

                align_area1_r=imshearotate(align_area1(:,:,aa),ann_last+ann);
                align_area2_r=imshearotate(align_area2(:,:,aa),-(ann_last+ann));
                %align_area1_r=align_area1_r(1+5:end-5,1+5:end-5,:);
                %align_area2_r=align_area2_r(1+5:end-5,1+5:end-5,:);
                for ii=scann
                    for jj=scann
                        align_area1_temp=real(s_ifft2(s_fft2(align_area1_r).*exp(-1i.*2.*pi.*((ii_last+ii)*coo1+(jj_last+jj)*coo2))));
                        align_area2_temp=real(s_ifft2(s_fft2(align_area2_r).*exp(+1i.*2.*pi.*((ii_last+ii)*coo1+(jj_last+jj)*coo2))));
                        %cost=sum(filter_border.*abs(align_area1_temp+align_area2_temp),'all');
                        cost=sum(abs(align_area1_temp+align_area2_temp),'all');
                        %no_neg=align_area1_temp+align_area2_temp;
                        %no_neg(no_neg<0)=0;
                        %cost=sum(filter_border.*abs(no_neg),'all');
                        if cost<cost_min
                            ii_opt=ii;
                            jj_opt=jj;
                            ann_opt=ann;
                            cost_min=cost;
                        end
                    end
                end
            end
            ii_last=ii_last+ii_opt;
            jj_last=jj_last+jj_opt;
            ann_last=ann_last+ann_opt;
            %error('stop')
            align_area1_out(:,:,aa)=real(s_ifft2(s_fft2(imshearotate(align_area1_out(:,:,aa),ann_opt)).*exp(-1i.*2.*pi.*(ii_opt*coo1+jj_opt*coo2))));
            align_area2_out(:,:,aa)=real(s_ifft2(s_fft2(imshearotate(align_area2_out(:,:,aa),-ann_opt)).*exp(+1i.*2.*pi.*(ii_opt*coo1+jj_opt*coo2))));
        end

        rawSino_rec_out(:,:,thetaInd1(aa))=real(s_ifft2(s_fft2(imshearotate(rawSino_rec_out(:,:,thetaInd1(aa)),ann_last)).*exp(-1i.*2.*pi.*(ii_last*coo1_big+jj_last*coo2_big))));
        rawSino_rec_out(:,:,thetaInd2(aa))=real(s_ifft2(s_fft2(imshearotate(rawSino_rec_out(:,:,thetaInd2(aa)),ann_last)).*exp(+1i.*2.*pi.*(ii_last*coo1_big-jj_last*coo2_big))));

        %raw_BG_rec_out(:,:,1)=real(s_ifft2(s_fft2(raw_BG_rec_out(:,:,1)).*exp(-1i.*2.*pi.*(ii_last*coo1_big+jj_last*coo2_big))));

        after=align_area1_out(:,:,aa)+align_area2_out(:,:,aa);
        display([num2str(aa) '/' num2str(size(align_area1,3)) ' ---> optimal : ' num2str(ii_last) '/' num2str(jj_last) '/' num2str(ann_last)]);
        imagesc(cat(2,before,after)); axis image; drawnow;
    end

    %%
    diff_before=rawSino_rec(:,:,1:end/2)+flip(rawSino_rec(:,:,end/2+1:end),2);
    diff_after=rawSino_rec_out(:,:,1:end/2)+flip(rawSino_rec_out(:,:,end/2+1:end),2);
    figure; sliceViewer(cat(2,diff_before,diff_after));
    %%
    clear align_area1 align_area2 align_area1_out align_area2_out;
    %%
    recon_aligned=true;

    crp_final=5;
    rawSino_rec_out_final=rawSino_rec_out(1+crp_final:end-crp_final,1+crp_final:end-crp_final,:);
    raw_BG_rec_out_final=raw_BG_rec_out(1+crp_final:end-crp_final,1+crp_final:end-crp_final,:);
    num_itter = 15;
    fieldReconInputParams = {num_itter, p, 'figShow',true,'DPCdata',struct([])};%DPCdata};
    if recon_aligned
        [yy, xx, Nang0] = size(rawSino_rec_out_final);
    else
        [yy, xx, Nang0] = size(rawSino_rec);
    end
    RI = zeros(yy,xx,xx,'single');
    %rawSino_rec_out_final=padarray(rawSino_rec_out_final,[0 round(size(rawSino_rec_out_final,2)/2)  0]);

    %rawSino_rec_out_final=cat(2,flip(rawSino_rec_out_final,2),rawSino_rec_out_final);

    %d1=((1:size(rawSino_rec_out_final,1))-floor(size(rawSino_rec_out_final,1)/2)+1)';d1=d1./max(d1(:));
    %d2=((1:size(rawSino_rec_out_final,2))-floor(size(rawSino_rec_out_final,2)/2)+1);d2=d2./max(d2(:));
    %rawSino_rec_out_final=real(s_fft2(s_ifft2(rawSino_rec_out_final)./(0.5+exp(-2000*(d1.^2+d2.^2)))));

    bytesPerSlice = (xx*Nang0) * 8; % 4 bytes for single
    gpu=gpuDevice;
    maxChunkSize = 1000;%floor( gpu.AvailableMemory / bytesPerSlice / 30 ); % last factor is a heuristic safety value.
    Npart = ceil( yy / maxChunkSize);
    chunks = round( linspace(1, yy+1, Npart+1) );

    for ii = 1:Npart
        fprintf( '[%s] chunk start: %03d/%03d\n', datetime,ii,Npart )
        thisChunk = chunks(ii) : chunks(ii+1)-1;
        if recon_aligned
            nmmm=1;
            validAngInd=[(1:nmmm:floor(size(rawSino_rec_out_final,3)/2)) (floor(size(rawSino_rec_out_final,3)/2)+1):nmmm:size(rawSino_rec_out_final,3)];
            RI(thisChunk,:,:) = sino2tomo(gpuArray(rawSino_rec_out_final(thisChunk,:,validAngInd)), angleLib(validAngInd), 0, bgWingSize, fieldReconInputParams, 0,raw_BG_rec_out_final(thisChunk,:));
            %RI(thisChunk,:,:) = sino2tomo(gpuArray(cat(3,circshift(rawSino_rec_out_final(thisChunk,:,validAngInd(1:end/2)),[0 0 0]),rawSino_rec_out_final(thisChunk,:,validAngInd(end/2+1:end)))), angleLib(validAngInd), 0, bgWingSize, fieldReconInputParams, 0,raw_BG_rec_out_final(thisChunk,:));
        else
            RI(thisChunk,:,:) = sino2tomo(gpuArray(rawSino_rec(thisChunk,:,validAngInd)), angleLib(validAngInd), 0, bgWingSize, fieldReconInputParams, 0,raw_BG_rec_out(thisChunk,:));
        end
        %RI(thisChunk,:,:) = sino2tomo(gpuArray(rawSino(thisChunk,:,validAngInd)), angleLib(validAngInd), tH.axisOffset, bgWingSize, fieldReconInputParams, 0);
        % pause;
    end

    figure, contrastVis(real(RI)); axis image; colorbar; %title(titleStr, 'Interpreter', 'none')
    figure, contrastVis(imag(RI)); axis image; colorbar;
%%
    save([MAIN_PATH name '_RI_result.mat'],'RI');

end