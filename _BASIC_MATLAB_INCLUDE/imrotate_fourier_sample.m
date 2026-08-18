function img=imrotate_fourier_sample(img,angle)

padd_size=round(size(img)/2);

img=fftshift(fft2(ifftshift(img)));
img=padarray(img,[padd_size(1) padd_size(2)],0,'both');
img=fftshift(ifft2(ifftshift(img)));
img=imrotate(img,angle);
img=fftshift(fft2(ifftshift(img)));
img=img(1+padd_size(1):end-padd_size(1),1+padd_size(2):end-padd_size(2),:,:);
img=fftshift(ifft2(ifftshift(img)));

end