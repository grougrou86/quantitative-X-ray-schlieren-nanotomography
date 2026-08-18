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

%% curve in the inner vs outer part of sample

inner=mean(gather(imag(registered_RI.*(rsize_filter<0.2))),[1 2 3]);
inner=squeeze(inner-mean(inner(1)));
shell=mean(gather(imag(registered_RI.*(rsize_filter<0.26).*(rsize_filter>0.2))),[1 2 3]);
shell=squeeze(shell-mean(shell(1)));
figure; plot(spectrum(:),shell(:));hold on; plot(spectrum(:),inner(:));

[max_shell,shell_fine,smooth]=get_max_graph(spectrum,shell);
[max_inner,inner_fine,smooth]=get_max_graph(spectrum,inner);

figure; plot(smooth,shell_fine(:));hold on; plot(smooth,inner_fine(:));

display(['Shell peak : ' num2str(max_shell)]);
display(['Core peak : ' num2str(max_inner)]);

%% Attempt 3D xanes
resize_for_xanes=0.1;

xanes_array=[];
for sp=1:size(registered_RI,4);
   cur_im= imresize3(imag(registered_RI(:,:,:,sp)),resize_for_xanes,'linear');
   if sp==1
       xanes_array=zeros([size(cur_im) size(registered_RI,4)],'single','gpuArray');
   end
   xanes_array(:,:,:,sp)=cur_im;
end
xanes_array=xanes_array-xanes_array(:,:,:,1);
max_pos_array=xanes_array(:,:,:,1).*0;

xanes_array_flat=reshape(xanes_array,[],size(xanes_array,4));
spectrum=gpuArray(single(spectrum));

xanes_array_flat=gather(xanes_array_flat);spectrum=gather(spectrum);max_pos_array=gather(max_pos_array);

pos_red=77593;%75049;
spectrum_red=[];spectrum_smooth_red=[];coo_smooth_red=[];
pos_green=103033;
spectrum_green=[];spectrum_smooth_green=[];coo_smooth_green=[];

for ii=1:length(max_pos_array(:))
    if mod(ii,1000)==0
        display([num2str(ii) ' / ' num2str(length(max_pos_array(:)))]);
    end
    [max_pos_array(ii),~,~]=get_max_graph(spectrum,squeeze(xanes_array_flat(ii,:)));
    if ii==pos_red 
        coo_red=spectrum;
        spectrum_red=squeeze(xanes_array_flat(ii,:));
        [~,spectrum_smooth_red,coo_smooth_red]=get_max_graph(spectrum,squeeze(xanes_array_flat(ii,:)));
    end
    if  ii==pos_green
        coo_green=spectrum;
        spectrum_green=squeeze(xanes_array_flat(ii,:));
        [~,spectrum_smooth_green,coo_smooth_green]=get_max_graph(spectrum,squeeze(xanes_array_flat(ii,:)));
    end
end


%{
 figure; orthosliceViewer(gather(xanes_array(:,:,:,4)));caxis([0 5e-7]);
 figure; orthosliceViewer(gather(max_pos_array));caxis([8371 8373]);
 figure; orthosliceViewer(gather(imresize3(max_pos_array,1/resize_for_xanes,'nearest')));caxis([8371 8373]);
%}
%%
id=1:numel(max_pos_array(:,:,:,1));
id=reshape(id,size(max_pos_array,1:3));
img_peak=imresize3(max_pos_array,size(registered_RI,1:3),'nearest');
id=imresize3(id,size(registered_RI,1:3),'nearest');



%img_peak=(img_peak-8371)/2;

fig=figure; orthosliceViewer(img_peak,Parent=fig);caxis([8371.25 8372.75]); colorbar;
fig=figure; orthosliceViewer(id,Parent=fig);

img_peak=(img_peak-8371.25)/1.5;
img_peak=1-img_peak;
img_peak=uint8((img_peak.*255)+1);

load('G:\20240911_BMOL\MATLAB\map_peak.mat','cmap');
figure; imagesc((255:-1:0)'./255.*ind2rgb(0:255,cmap)); 
map=cmap;%[256:-1:1;1:256;0.*(1:256)]'./256;
%map=map./max(map,[],2);
img_peak=reshape(ind2rgb(img_peak(:,:),map),[size(img_peak) 3]);



disp_intensity=imag(mean(registered_RI(:,:,:,:),4));

disp_intensity=disp_intensity./1e-6;

disp_intensity(disp_intensity<0)=0;
disp_intensity(disp_intensity>1)=1;

fig=figure; orthosliceViewer(disp_intensity,Parent=fig);
fig=figure; orthosliceViewer(disp_intensity.*img_peak,Parent=fig);
%%

figure; hold on; 

shift=0;

min_r=min(spectrum_smooth_red(:));
max_r=max(spectrum_smooth_red(:));
plot(coo_smooth_red,(spectrum_smooth_red-min_r)./(max_r-min_r)-shift,'-','Color','red'); 
plot(coo_red,(spectrum_red-min_r)./(max_r-min_r)-shift,'o','Color','red','MarkerSize',7); 
ylim([-0.2 1.2]);

min_g=min(spectrum_smooth_green(:));
max_g=max(spectrum_smooth_green(:));
plot(coo_smooth_green,(spectrum_smooth_green-min_g)./(max_g-min_g)+shift,'-','Color','green');
plot(coo_green,(spectrum_green-min_g)./(max_g-min_g)+shift,'o','Color','green','MarkerSize',7);
