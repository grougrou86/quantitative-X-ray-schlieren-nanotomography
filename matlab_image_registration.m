img1=single(imread('F:\2024\20240911_BMOL\angle_test\TOMO_00000.tif'));
img2=single(imread('F:\2024\20240911_BMOL\angle_test\TOMO_00450.tif'));
bg=single(imread('F:\2024\20240911_BMOL\angle_test\FF_00000.tif'));

%img1=single(imread('F:\2023\20230728 LKRHerve\tomodata\gradient_nmc1\0deg_1s.tif'));
%img2=single(imread('F:\2023\20230728 LKRHerve\tomodata\gradient_nmc1\180deg_1s.tif'));
%bg=single(imread('F:\2023\20231224 BMOL\tomodata\gradient_nmc1_1s_3_proj\gradient_nmc1_1s_3_ff_0001.tif'));

%img1=single(imread('F:\2023\20231224 BMOL\tomodata\gradient_nmc1_1s_3_proj\gradient_nmc1_1s_3_0001.tif'));
%img2=single(imread('F:\2023\20231224 BMOL\tomodata\gradient_nmc1_1s_3_proj\gradient_nmc1_1s_3_0901.tif'));
%bg=single(imread('F:\2023\20231224 BMOL\tomodata\gradient_nmc1_1s_3_proj\gradient_nmc1_1s_3_ff_0001.tif'));

save_name='F:\2023\20231224 BMOL\tomodata\temp_for _computation\register_field';


% bg substraction

%min_inv=mean(bg(:))/10;
%bg(bg<min_inv)=min_inv;


img1=img1./mean(img1,2);
img2=img2./mean(img2,2);
bg=bg./mean(bg,2);
%%
img1=img1./bg;
img2=img2./bg;

img2=flip(img2,2);


figure; imagesc(img1-img2);

img1_0=img1;
img2_0=img2;

half_rotation=0;

%%
for ii=1:2

    img1=imrotate(img1_0,-half_rotation,'crop');
    img2=imrotate(img2_0,half_rotation,'crop');

    img1=mcrop(img1,size(img1)-[100 100]);
    img2=mcrop(img2,size(img2)-[100 100]);

    crop=round(size(img1,1)/50);

    high_pass_1=s_fft2(img1);
    high_pass_1(floor(end/2)+(-crop:crop),floor(end/2)+(-crop:crop))=0;
    high_pass_1=real(s_ifft2(high_pass_1));
    high_pass_1=mcrop(high_pass_1,size(high_pass_1)-[100 100]);

    high_pass_2=s_fft2(img2);
    high_pass_2(floor(end/2)+(-crop:crop),floor(end/2)+(-crop:crop))=0;
    high_pass_2=real(s_ifft2(high_pass_2));
    high_pass_2=mcrop(high_pass_2,size(high_pass_2)-[100 100]);
    %% register
    [optimizer,metric] = imregconfig("multimodal");
    [registeredDefault,R_reg,tform] = imregister2(high_pass_1,-high_pass_2,"rigid",optimizer,metric);
    tform
    if ii==1
        half_rotation=tform.RotationAngle/2;
    end
    %[registeredDefault,R_reg,tform] = imregister2(high_pass_1,-high_pass_2,"translation",optimizer,metric);
    %tform = imregcorr(high_pass_1,-high_pass_2,"rigid");
    %moving = imwarp(registeredDefault,tform);

    %figure; imshowpair(moving,high_pass_2,'montage');
    figure; imshowpair(registeredDefault,high_pass_2,'montage');
    %%
    Rfixed = imref2d(size(img2));
    Rmoving = imref2d(size(img1));
    [img1Reg,Rreg] = imwarp(img1,Rmoving,tform,'OutputView',Rfixed);

    gradient=mcrop(img2-img1Reg,[800 800]);

    figure; imagesc(gradient);

    diff_img_2=mcrop(img2,size(img2)-[200 200]);
    diff_img_1=mcrop(img1Reg,size(img1Reg)-[200 200]);
    figure; imagesc(diff_img_1-diff_img_2);axis image; colormap gray;

end

%save(save_name,'diff_img_1','diff_img_2','tform','half_rotation')
