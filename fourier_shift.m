function Array_out=fourier_shift(Array_in,d1,d2,d3)
if strcmp(class(Array_in),'gpuArray')
    [coor2,coor1,coor3] = meshgrid(gpuArray(single(1:size(Array_in,2))),gpuArray(single(1:size(Array_in,1))),gpuArray(single(1:size(Array_in,3))));
else
    [coor2,coor1,coor3] = meshgrid((single(1:size(Array_in,2))),(single(1:size(Array_in,1))),(single(1:size(Array_in,3))));
end
centre1=floor(size(coor1,1)/2)+1;
centre2=floor(size(coor2,2)/2)+1;
centre3=floor(size(coor3,3)/2)+1;
coor1=coor1-centre1;
coor2=coor2-centre2;
coor3=coor3-centre3;

Array_out=fftshift(fftn(ifftshift(Array_in)));
filter=exp(-1i*(2*pi)*(d1.*coor1./size(Array_in,1)+d2.*coor2./size(Array_in,2)+d3.*coor3./size(Array_in,3)));
Array_out=Array_out.*filter;
Array_out=fftshift(ifftn(ifftshift(Array_out)));
end