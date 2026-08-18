function [Eout, DPCdata] =field_recon_routine_v21(sino1, sino2,BG, num_itter, pr, varargin)

symmetrize=false;true;
if symmetrize
    sino1=cat(2,sino1,sino2);
    sino2=cat(2,sino2,sino1(:,1:end/2,:));
    BG=cat(2,BG,BG);
    %sino1=cat(2,sino1,flip(sino1,2));
    %sino2=cat(2,sino2,flip(sino2,2));
    %BG=cat(2,BG,flip(BG,2));
end

BG=permute(BG,[2 1]);
sino1=permute(sino1,[2 1 3]);
sino2=permute(sino2,[2 1 3]);

p = inputParser();
addParameter(p, 'figShow', true, @islogical);
addParameter(p, 'wienerConst', 1e-8, @isscalar);%1e-8
addParameter(p, 'deconvConst', 0.8, @isscalar);
addParameter(p, 'DPCdata', struct([]), @isstruct);
addParameter(p, 'aberration', 1);

% parse
parse(p, varargin{:});
figShow = p.Results.figShow;
wienerConst = p.Results.wienerConst;
deconvConst = p.Results.deconvConst;
DPCdata = p.Results.DPCdata;
aberation=p.Results.aberration;

%% asssertions
assert( all(size(sino1) == size(sino2)) )

if isempty(DPCdata)
    xx = size(sino1, 1);
    assert( isfield(pr, 'ZP_focal_length') );


    %% basic parameters
    wl = pr.wavelength;
    imgPix = pr.pixel_size;
    % sample_lens_dist = params_recon.sample_lens_dist;
    FimgPix = pr.ZP_focal_length .* wl / (xx.*imgPix); %size of on pixel in the fourier space in 'um' as measured at the beam block

    %% make the pattern and pupil

    %%% cutoff pupil
    pupil   = zeros(xx*4,1,1,2,'single'); % [yy,xx, ]
    %mx = floor(xx/2)+1;
    mx = floor(xx*2)+1;
    crp=0;%floor(xx/3);
    pupil(1+crp:mx,:,1,1) = 1;
    pupil(mx:end-crp,:,1,2) = 1;
    pupil=s_fft(pupil,1);
    pupil=pupil(floor(end/2)-floor(end/4/2)+(1:end/4),:,:,:)/4;
    pupil=s_ifft(pupil,1);
    pupil=pupil;%.*aberation;


    mx = floor(xx/2)+1;

    %%% measured illumination pattern
    % plot(p.fourier_position, p.fourier_intensity)
    % pattern = zeros(1,xx,'single');

    pattern   = zeros(xx,1,1,2,'single'); % [yy,xx, ]

    % figure,plot(kx0, pr.fourier_intensity)

    if isfield(pr,'gaussFit')% && false
        % kx0 = (pr.fourier_position)*1e-3;
        kx1 = ((1:xx)'-mx)*FimgPix; % 0 --> cutoff position
        boundary  = (pr.gaussFit.b + 3*pr.gaussFit.c *[-1,1]); % 3 sigma (99.7%)

        xvec1 = (+kx1 + pr.cutoff)*1e-3;
        xvec2 = (-kx1 + pr.cutoff)*1e-3;
        pattern(1,:,1,1) = pr.gaussFit(xvec1);
        pattern(1,:,1,2) = pr.gaussFit(xvec2);
        pattern = pattern -  pr.gaussFit.d; % offset zeroing

        validInd1 = xvec1 < boundary(2) & xvec1 > boundary(1);
        validInd2 = xvec2 < boundary(2) & xvec2 > boundary(1);

        pattern(~validInd1,1,1,1) = 0;
        pattern(~validInd2,1,1,2) = 0;
        % figure, plot(kx1,pattern(:,:,:,1),kx1,pattern(:,:,:,2));

    else
        kx0 = pr.fourier_position - pr.cutoff;%+1;
        %pr.fourier_position %-1711-1613
        %pr.cutoff%-1662
        %error('remove offset')
        kx1 = ((1:xx)'-mx)*FimgPix;

        pattern(:,:,1,1) = interp1( kx0, pr.fourier_intensity, kx1, 'linear',0 );
        pattern(:,:,1,2) = interp1( kx0, pr.fourier_intensity, -kx1, 'linear',0 );
        pattern(pattern < 0) = 0;
    end
    % figure, plot(kx1,pattern(:,:,:,1),kx1,pattern(:,:,:,2));


    %% compute the PSF
    %pupil=imgaussfilt(pupil,100);
    filtered_pattern = pattern.*pupil;

    if figShow
        figure('Name','pattern and filtered_pattern'); plot(pattern(:,:,:,1));hold on; plot((filtered_pattern(:,:,:,1).*0.95));
    end

    filtered_pattern = filtered_pattern./(sum(abs(filtered_pattern), [1 2]));
    % plot(kx1,filtered_pattern(:,:,:,1),kx1,filtered_pattern(:,:,:,2));

    PSF = conj(s_ifft(pupil,1)) .* s_ifft(filtered_pattern,1) .* prod(size(pattern,[1,2])); %
    % plot(kx1, abs(s_fft(PSF(:,:,:,1),2)), kx1, abs(s_fft(PSF(:,:,:,2),2)));

    H=2*real(s_fft(real(PSF),1)); % diffraction psf
    P=2*imag(s_fft(imag(PSF),1)); % diffraction psf
    

    % plot(kx1, H(:,:,:,1), kx1, H(:,:,:,2));
    % plot(kx1, P(:,:,:,1), kx1, P(:,:,:,2));
    %size(H)
    %% to plot the psf etc..
    %figure; plot(((1:xx)'-mx)*1/(xx*imgPix),H(:,1)); hold on; plot(((1:xx)'-mx)*1/(xx*imgPix),P(:,1)); xlim([-11.4 11.4])
    %figure; plot(((1:xx)'-mx)*1/(xx*imgPix),pupil(:,1));hold on; plot(((1:xx)'-mx)*1/(xx*imgPix),filtered_pattern(:,1));xlim([-11.3 11.3])
    %error('stop')
    %% inversion
    deconv_range = (H(:,:,1) > deconvConst);

    %wienerConst=1e-3;

    Inv_1=sum(H(:,:,:).*conj(H(:,:,:))+wienerConst,3);
    Inv_2=sum(H(:,:,:).*conj(P(:,:,:)),3);
    Inv_3=sum(P(:,:,:).*conj(H(:,:,:)),3);
    Inv_4=sum(P(:,:,:).*conj(P(:,:,:))+wienerConst,3);

    det=1./(Inv_1.*Inv_4-Inv_2.*Inv_3);

    invert_mat  = zeros([size(H,1) size(H,2) 2 2 ],'single');
    invert_mat2 = zeros([size(H,1) size(H,2) 2 2 ],'single');

    invert_mat(:,:,1,1)=conj(H(:,:,1).*det.*deconv_range);
    invert_mat(:,:,1,2)=conj(H(:,:,2).*det.*deconv_range);
    invert_mat(:,:,2,1)=conj(P(:,:,1).*det.*deconv_range);
    invert_mat(:,:,2,2)=conj(P(:,:,2).*det.*deconv_range);

    invert_mat2(:,:,1,1)=conj(+Inv_4);
    invert_mat2(:,:,1,2)=conj(-Inv_2);
    invert_mat2(:,:,2,1)=conj(-Inv_3);
    invert_mat2(:,:,2,2)=conj(+Inv_1);

    invert_mat=page_time_end(invert_mat2,invert_mat);
    invert_mat=sum(invert_mat,3);
    invert_mat=single(reshape(invert_mat,size(invert_mat,1),size(invert_mat,2),1,[]));
    % plot(kx1, invert_mat(:,:,:,1), kx1, invert_mat(:,:,:,2));


    %% dataOut
    DPCdata = struct;
    DPCdata.invert_mat = invert_mat;
    DPCdata.pattern = pattern;
    DPCdata.pupil = pupil;
else

    invert_mat = DPCdata.invert_mat;
    pattern    = DPCdata.pattern;
    pupil      = DPCdata.pupil;
end


%% DPCrecon iteration
%inputData = gpuArray(single(cat(4,sino2,sino1)));
inputData = sino1.*reshape([0 1],1,1,1,[])+sino2.*reshape([1 0],1,1,1,[]);
if ~isa(sino1,'ADNode')
    clear sino1 sino2
end

BG  = ifftshift(BG,1);
inputData  = ifftshift(inputData,1);
invert_mat = ifftshift(gpuArray(invert_mat),1);
pattern    = ifftshift(gpuArray(pattern),1);
pupil      = ifftshift(gpuArray(pupil),1);

%[yy, xx, Nang180] = size(inputData,[1,2,3]);
sz_input_data= size(inputData);
yy=sz_input_data(1);
xx=sz_input_data(2);
Nang180=sz_input_data(3);
Eout = zeros(yy, xx, Nang180, 'single', 'gpuArray');
sim_intensity_image = ones(yy, xx, Nang180, 2, 'single', 'gpuArray');
erro = nan(1, num_itter);
% bytesPerSlice = (xx*Nang180) * 8; % 8 bytes for complex single
% gpu=gpuDevice;
% maxChunkSize = floor( gpu.AvailableMemory / bytesPerSlice / 10 ); % last factor is a heuristic safety value.
% Npart = ceil( yy / maxChunkSize);
% chunks = round( linspace(1, yy+1, Npart+1) );
if figShow
    figure('Name','Reconstructed field');
end

err_pos=[];

pseudo_meth=3;%pseudo newton J_0 = 1 / pseudo newton J_0 + intensity correction = 3 / Newton = 4

inputData=inputData-mean(inputData,[1])+mean(inputData,[1,4]);
errorPoints=false;
for kk = 1:num_itter
    tic;
    %res = ifft( sum(invert_mat.*fft(inputData - sim_intensity_image,[],1),4),[],1);
    switch pseudo_meth
        case 1
            res = ifft( sum(fft(inputData - sim_intensity_image,[],1).*invert_mat,4),[],1);
        case 2
            error('use 3 not 2. 2 in unstable')
            res = ifft( sum(fft((inputData - sim_intensity_image)./mean(sim_intensity_image,4),[],1).*invert_mat,4),[],1);
            %res = ifft( sum(fft((inputData - sim_intensity_image).*abs(exp(-2*real(Eout))),[],1).*invert_mat,4),[],1);
        case 3
            %res = abs(exp(-2*real(Eout))).*ifft( sum(fft((inputData - sim_intensity_image),[],1).*invert_mat,4),[],1);
            res = 1./max(sim_intensity_image,[],4).*ifft( sum(fft((inputData - sim_intensity_image),[],1).*invert_mat,4),[],1);

            %res = imgaussfilt(abs(exp(-2*real(Eout))),4).*ifft( sum(fft((inputData - sim_intensity_image),[],1).*invert_mat,4),[],1);
        otherwise
            error=(inputData - sim_intensity_image);
            ident=diag(ones(2*size(error,1),1,'single','gpuArray'));
            %alpha_filter=diag(cat(1,ones(size(error,1),1,'single','gpuArray'),zeros(size(error,1),1,'single','gpuArray')));
            %phi_filter=diag(cat(1,zeros(size(error,1),1,'single','gpuArray'),ones(size(error,1),1,'single','gpuArray')));
            alpha_filter=diag(ones(size(error,1),1,'single','gpuArray'));
            alpha_filter=cat(1,alpha_filter,alpha_filter);

            phi_filter  =cat(2,0.*alpha_filter,1.*alpha_filter);
            alpha_filter=cat(2,1.*alpha_filter,0.*alpha_filter);


            fft_mat=fft(diag(ones(size(error,1),1,'single','gpuArray')),[],1);
            fft_mat=cat(2,cat(1,fft_mat,0.*fft_mat),cat(1,0.*fft_mat,fft_mat));
            ifft_mat=ifft(diag(ones(size(error,1),1,'single','gpuArray')),[],1);
            ifft_mat=cat(2,cat(1,ifft_mat,0.*ifft_mat),cat(1,0.*ifft_mat,ifft_mat));


            deconv_range_filter=diag(gpuArray(single(cat(1,ifftshift(deconv_range),ifftshift(deconv_range)))));
            deconv_range_filter=pagemtimes(ifft_mat,pagemtimes(deconv_range_filter,fft_mat));
            %filter_def
            n=size(ifft_mat,1);

            subdiv=4;
            res=0.*Eout;
            for qq=1:size(Eout,3)
                display(['qq = ' num2str(qq)]);
                for ii=1:subdiv
                    display(['ii = ' num2str(ii)]);
                    chunk=round((ii-1)*size(error,2)/subdiv)+1:round(ii*size(error,2)/subdiv);
                    p=length(chunk);
                    %p=size(Eout,3);
                    diag_field=zeros(n^2,p,'single','gpuArray');
                    diag_patt_exp=zeros(n^2,1,'single','gpuArray');
                    temp_B_diag=zeros(n^2,p,'single','gpuArray');
                    diag_filter=diag(pupil(:));

                    x=error(:,chunk,qq,:);
                    %x=permute(x,[ 1 4 3 2]);
                    x=permute(x,[ 1 4 2 3]);
                    x=reshape(x,size(x,1)*size(x,2),1,size(x,3),size(x,4));
                    %field def
                    field=(cat(1,exp(Eout(:,chunk,qq)),exp(Eout(:,chunk,qq))));field=reshape(field,size(field,1),1,[]);

                    diag_field=reshape(diag_field,n^2,p);
                    diag_field(1:n+1:end,:)=reshape(field,n,p);
                    diag_field=reshape(diag_field,n,n,p);

                    pattern_norm = sqrt(pattern./sum(abs(pupil.*pattern), [1,2])) .*size(pattern,1);
                    %pattern_kk = pattern_norm(:,:,:,kk);
                    nonzero_ill_1 = find( pattern(:,:,:,1)>0 );
                    nonzero_ill_2 = find( pattern(:,:,:,1)>0 );
                    Nill = max(length(nonzero_ill_1),length(nonzero_ill_2));
                    A=0;

                    for pp = 1:Nill
                        %Nill
                        %pp
                        temp = zeros(size(Eout,1),1,1,2,'single','gpuArray');
                        if(pp<=length(nonzero_ill_1))
                            temp(nonzero_ill_1(pp)) = pattern_norm(nonzero_ill_1(pp)); % 1d
                        end
                        if(pp<=length(nonzero_ill_2))
                            temp(size(Eout,1)+nonzero_ill_2(pp)) = pattern_norm(size(Eout,1)+nonzero_ill_2(pp)); % 1d
                        end
                        %temp = ifft(temp,[],1);
                        temp=temp(:);
                        temp=pagemtimes(ifft_mat,temp);

                        diag_patt_exp=reshape(diag_patt_exp,n^2,1);
                        diag_patt_exp(1:n+1:end,:)=reshape(temp,n,1);
                        diag_patt_exp=reshape(diag_patt_exp,n,n,1);
                        temp_A=pagemtimes(ifft_mat,pagemtimes(diag_filter,pagemtimes(fft_mat,pagemtimes(diag_patt_exp,diag_field))));
                        %temp_B=pagemtimes(ifft_mat,pagemtimes(diag_filter,pagemtimes(fft_mat,pagemtimes(diag_patt_exp,field))));
                        temp_B=pagemtimes(temp_A,field*0+1);

                        temp_B_diag=reshape(temp_B_diag,n^2,p);
                        temp_B_diag(1:n+1:end,:)=reshape(temp_B,n,p);
                        temp_B_diag=reshape(temp_B_diag,n,n,p);
                        temp_A=pagemtimes(conj(temp_B_diag),temp_A);
                        A=A+temp_A;
                    end

                    A_final=A;
                    A_final=(pagemtimes(A_final,1.*alpha_filter+1i.*phi_filter)+pagemtimes(conj(A_final),1.*alpha_filter-1i.*phi_filter));

                    conj_val=conj(pagefun(@transpose,A_final));
                    A_inv=pagemtimes(pagefun(@inv,pagemtimes(conj_val,A_final)+ident*wienerConst*1000),conj_val);
                    A_inv=pagemtimes(A_inv,deconv_range_filter);
                    y=pagemtimes(A_inv,x);
                    y=y(1:end/2,:)+1i.*y(end/2+1:end,:);
                    res(:,chunk,qq)=y;

                end
            end
            %{
            res = ifft( sum(fft(inputData - sim_intensity_image,[],1).*invert_mat,4),[],1);
            res_cat=cat(1,real(res),imag(res));
            res_cat=res_cat(:,chunk,1);
            test_for=pagemtimes(A_final,res_cat);
            %}
    end
    Eout_old=Eout;
    Eout = 1*res + Eout;%0.75 has slightly better stability
    %Eout(real(Eout)>0)=Eout(real(Eout)>0)-real(Eout(real(Eout)>0));

    if kk<num_itter

        sim_intensity_image = DPC_simulation( exp(Eout), pattern, pupil ,BG);
        
        errorPoints_last=errorPoints;
        errorPoints = (abs(sim_intensity_image - inputData) > 10)|isnan(abs(sim_intensity_image - inputData));
        errorPoints=errorPoints|errorPoints_last;
        if any(errorPoints,'all')

            errorInds = find(errorPoints(:));
            [~, ~, errorAngle, ~] = ind2sub(size(sim_intensity_image) , errorInds);
            %errorAngle
            display(['RMS to high for fields : ' num2str(gather(unique(errorAngle(:)')))]);
            err_pos=unique(errorAngle(:));
            %errorAngStr = sprintf('%05d, ', unique(errorAngle));
            %error(' RMSE is too high.. errorAngles: \n%s', errorAngStr)
        end


    end

    if ~isempty(err_pos)
        Eout(:,:,err_pos,:)=0;
        sim_intensity_image(:,:,err_pos,:)=0;
    end

    if kk<num_itter
        rmse_temp=rmse(sim_intensity_image, inputData, [1 2 4]);
        rmse_temp(:,:,err_pos,:)=0;
        curr_erro=mean(rmse_temp,'all');
        if kk>1 && erro(kk-1)<curr_erro
            Eout=Eout_old;
            break;
        end
        erro(kk) = curr_erro;
        
    end

    if figShow
        nn = 2; mm = 2;
        cropp_out=fftshift(Eout(:,:,1),1);
        if symmetrize
            cropp_out=cropp_out(1:floor(end/2),:,:);
        end
        cropp_out=cropp_out-1i.*imag(mean(cropp_out(1:floor(end/100),:,:),1));
        ax=subplot(nn,mm,1); imagesc(real(cropp_out)); axis image; colorbar; title('log(A)'); colormap(ax,gray);
        ax=subplot(nn,mm,2); imagesc(imag(cropp_out)); axis image; colorbar; title('phase');colormap(ax,turbo);
        ax=subplot(nn,mm,3); imagesc(fftshift(inputData(:,:,1,1)-sim_intensity_image(:,:,1,1), 1));  axis image; colorbar; title('diffmap')
        ax=subplot(nn,mm,4); semilogy(erro); axis square; ylabel('rmse'); xlabel('iteration')
        drawnow;
    end
    fprintf( '[%s] iter: %03d/%03d ... rmse: %.5f\n',datetime, kk,num_itter, erro(kk) )
    toc;
end

Eout = fftshift(Eout,1);
if symmetrize
    Eout=Eout(1:floor(end/2),:,:);
end
Eout=permute(Eout,[2 1 3]);
%stop 
end

function intensity_image = DPC_simulation(ground_truth, pattern, pupil,BG)

%%% ground_truth =
xx = size(pattern, 1);
sz=size(ground_truth);
if length(sz)<3
sz=[sz 1];
end
sz=sz(1:3);

assert(sz(1) == xx)


% pattern = pattern./mean(pupil.*pattern,[1, 2]);
% figure, plot(pattern(:,:,:,1));
%pupil=(pupil);
pattern = sqrt(pattern./sum(abs(pupil.*pattern), [1,2])) .*size(pattern,1);
intensity_image = zeros([ sz, size(pattern,4) ],'like', real(ground_truth(1)));
%size(pattern)
%size(BG)

for kk = 1:size(pattern,4)
    pattern_kk = pattern(:,:,:,kk);
    pupil_kk   = pupil(:,:,:,kk);
    nonzero_ill = find( pattern_kk>0 );
    Nill = length(nonzero_ill);

    for pp = 1:Nill

        temp = zeros(xx,1,'like',ground_truth);
        temp(nonzero_ill(pp)) = pattern_kk(nonzero_ill(pp)); % 1d
        temp = ifft(temp,[],1);
        temp = temp.*ground_truth.*sqrt(BG);%%BG norm
        temp = fft(temp,[],1);
        temp = temp.*pupil_kk;
        temp = ifft(temp,[],1);

        intensity_image(:,:,:,kk) = intensity_image(:,:,:,kk) + abs(temp).^2./BG;%%BG norm
    end

    % for pp = 1:Nill
    %
    %     temp = zeros(1,xx,'like',ground_truth);
    %     temp(nonzero_ill(pp)) = pattern_kk(nonzero_ill(pp)); % 1d
    %     temp = s_ifft(temp,2);
    %     temp = temp.*ground_truth;
    %     temp = s_fft(temp,2);
    %     temp = temp.*pupil_kk;
    %     temp = s_ifft(temp,2);
    %
    %     intensity_image(:,:,:,kk) = intensity_image(:,:,:,kk) + abs(temp).^2;
    % end
end
end

function out=s_fft(in,n)
out=fftshift(fft(ifftshift(in),[],n));
end

function out=s_ifft(in,n)
out=fftshift(ifft(ifftshift(in),[],n));
end
%figure; sliceViewer(gather(imag(tot_res)));
%figure; sliceViewer(gather(real(tot_res)));

