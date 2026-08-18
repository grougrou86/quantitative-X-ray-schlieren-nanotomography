scale=1.4;
angl=5;

img=imread("cameraman.tif");
img=gpuArray(single(padd_crop_to_fit(img,size(img).*2)));
figure; imagesc(img);axis image

img2=imshearotate(imresize_fft(img,scale),angl);


figure; imagesc(real(img2));axis image
