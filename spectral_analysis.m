MAIN_PATH='G:\20240911_BMOL\NCM2'
addpath('G:\20240911_BMOL\MATLAB\998 codes')
%MAIN_PATH='G:\20240911_BMOL\discharge_sample_2_with_filter'
%MAIN_PATH='G:\20240911_BMOL\charged';
files = dir([MAIN_PATH '\*result_new_reg.mat']);
name_list={files.name};
%%
spectrum=[];
RI_spec=[];
sino_spec=[];
kk=1;
name_list=name_list(5:1:end-1);
name_list=name_list([1:3 5:end]);

for name=name_list
    nums=str2double(extract(name{1}, digitsPattern));
    spectrum(kk)=nums(1);
    load([MAIN_PATH '\' name{1}]);
    [MAIN_PATH '\' name{1}]
    if isempty(RI_spec)
        RI_spec=complex(ones([size(RI) length(name_list)] ,'single'));
        %sino_spec=ones([size(rawSino_rec_out_final) length(name_list)] ,'single');
        
    end
    %sino_spec(:,:,:,kk)=rawSino_rec_out_final;
    RI_spec(:,:,:,kk)=RI;
    kk=kk+1;
    figure; imagesc(imag(gather(squeeze(gather(RI(:,:,floor(end/2)+1,:))))));
end

%% registraction 

registered_RI=RI_spec;



ref=RI_spec(:,:,:,1);
co=make_coo(size(ref));

filter_fourier=(1-exp(-500.*(co{3}.^2+co{2}.^2)));
img_filter=((co{3}.^2+co{2}.^2)<0.48.^2).*(abs(co{1})<0.49);

%ref=s_fftn(s_ifftn(s_fftn(gpuArray(single(real(ref)))).*filter_fourier).*img_filter);
ref=s_fftn(s_ifftn(s_fftn(gpuArray(single(real(ref)))).*filter_fourier).*img_filter);
%%


[optimizer, metric]  = imregconfig('monomodal');

scale_px=-10:1:10;

wait(gpuDevice());
for kk=1:size(RI_spec,4)
    kk
    curr_O=gpuArray(single((RI_spec(:,:,:,kk))));
    mx_val=0;
    for ss=scale_px
        curr=s_fftn(curr_O);
        curr=msize(curr,size(curr_O)+ss);
        curr=s_ifftn(curr);
        curr=msize(curr,size(curr_O));
        
        curr_scl=curr;
        
        curr=s_fftn(real(curr));
        curr=curr.*filter_fourier;
        curr=s_ifftn(curr);
        
        %curr_scl=curr;
        
        %curr=curr+10
        curr=curr.*img_filter;
        %error('stop')
        
        %curr=curr-mean(curr(:));
        curr=curr./sqrt(mean(curr(:).^2));
        
        curr=s_fftn(curr);
        coor=s_ifftn(curr.*conj(ref));
        [mxv,idx] = max(abs(coor(:)));
        
        if mxv>mx_val
            mx_val=mxv;
            [r,c,p] = ind2sub(size(coor),idx);
            ids=[r,c,p];
            ids=ids-floor(size(coor)/2)-1;
            ss_opt=ss;
            curr_opt=curr_scl;
        end
    end
    display(['result : scale = ' num2str(ss_opt) ' ; position = ' num2str(ids)])
    registered_RI(:,:,:,kk)=gather(circshift(curr_opt,-ids));
end
registered_RI_0=registered_RI;
%%
norm_RI= (registered_RI);
norm_RI=imag(norm_RI).*1i + real(norm_RI)./mean(real(norm_RI(end-50:end,:,:,:)),[1 2 3]).*mean(real(norm_RI(end-50:end,:,:,1)),[1 2 3]);
%%
registered_RI=s_fft3(registered_RI_0);
rsize_filter=sqrt(co{1}.^2+co{2}.^2+co{3}.^2);
registered_RI=registered_RI.*(exp(-2.355*(1*rsize_filter).^2));%FWHM at res/1
registered_RI=s_ifft3(registered_RI);
%%
figure;
yyaxis left;
mean_vals=squeeze(mean(registered_RI(round(2*end/6):round(4*end/6),round(2*end/6):round(4*end/6),round(2*end/6):round(4*end/6),:),[1 2 3]));
plot(spectrum(:),imag(mean_vals(:)),'-',Marker='o');
ylim([0 1.2e-6]);
hold on; 
yyaxis right;
plot(spectrum(:),real(mean_vals(:)),'-',Marker='o');
ylim([0 2e-5]);
%%
img=registered_RI(:,round(end/2),:,[1 3 4 6 8 10]);
img=img(:,:);
figure; 
subplot(2,1,1);
imagesc(real(img)); caxis([0 2e-5]); colormap(1-gray); axis image;
subplot(2,1,2);
imagesc(imag(img)); caxis([0 0.2e-5]); colormap(1-gray);axis image;
%%
figure; imagesc(1-2*cat(3,1.*(0:0.001:1)+0.*(0:0.001:1)',0.5.*(0:0.001:1)+0.5.*(0:0.001:1)',0.*(0:0.001:1)+1.*(0:0.001:1)'))
figure; imagesc(1-cat(3,1.*(0:0.001:1)+0.*(0:0.001:1)',0.5.*(0:0.001:1)+0.5.*(0:0.001:1)',0.*(0:0.001:1)+1.*(0:0.001:1)'))
%%
dd=interp1(spectrum,imag(mean_vals(:)),spectrum(1):1:spectrum(end));
padd=400;
dd=cat(2,linspace(dd(end),dd(1),padd),dd);
%dd=padarray(dd,[0 size(dd,2)],dd(1),'pre');

%dd=padarray(dd,[0 size(dd,2)],dd(end),'post');

dd=s_fft(dd,2);dd(1:floor(end/2)+1)=0;dd=2*s_ifft(dd,2);
dd=dd(padd:end);
figure; plot(real(dd)); hold on; plot(imag(dd));