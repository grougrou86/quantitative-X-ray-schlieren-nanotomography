function OTF=make_bright_field_OTF(size_img,resolution,NA,wavelength)


utility=struct;
%image space
utility.image_space=struct;
utility.image_space.res=cell(2);
mlt=4;
utility.image_space.res{1}=resolution(1)/mlt;
utility.image_space.res{2}=resolution(2)/mlt;
utility.image_space.size=cell(2);
utility.image_space.size{1}=size_img(1)*mlt;
utility.image_space.size{2}=size_img(2)*mlt;
utility.image_space.coor=cell(2);
utility.image_space.coor{1}=single((1:utility.image_space.size{1})-(floor(utility.image_space.size{1}/2)+1));
utility.image_space.coor{2}=single((1:utility.image_space.size{2})-(floor(utility.image_space.size{2}/2)+1));
utility.image_space.coor{1}=utility.image_space.coor{1}.*utility.image_space.res{1};
utility.image_space.coor{2}=utility.image_space.coor{2}.*utility.image_space.res{2};
utility.image_space.coor{1}=reshape(utility.image_space.coor{1},[],1,1);
utility.image_space.coor{2}=reshape(utility.image_space.coor{2},1,[],1);
%fourier space
utility.fourier_space=struct;
utility.fourier_space.res=cell(2);
utility.fourier_space.res{1}=1/(utility.image_space.res{1}*utility.image_space.size{1});
utility.fourier_space.res{2}=1/(utility.image_space.res{2}*utility.image_space.size{2});
utility.fourier_space.size=cell(2);
utility.fourier_space.size{1}=utility.image_space.size{1};
utility.fourier_space.size{2}=utility.image_space.size{2};
utility.fourier_space.coor=cell(2);
utility.fourier_space.coor{1}=single((1:utility.image_space.size{1})-(floor(utility.image_space.size{1}/2)+1));
utility.fourier_space.coor{2}=single((1:utility.image_space.size{2})-(floor(utility.image_space.size{2}/2)+1));
utility.fourier_space.coor{1}=utility.fourier_space.coor{1}.*utility.fourier_space.res{1};
utility.fourier_space.coor{2}=utility.fourier_space.coor{2}.*utility.fourier_space.res{2};
utility.fourier_space.coor{1}=reshape(utility.fourier_space.coor{1},[],1,1);
utility.fourier_space.coor{2}=reshape(utility.fourier_space.coor{2},1,[],1);
utility.fourier_space.coorxy=sqrt(...
    (utility.fourier_space.coor{1}).^2+...
    (utility.fourier_space.coor{2}).^2);
%other
utility.lambda=wavelength;
utility.k0=1/wavelength;
utility.kmax=NA/wavelength;
utility.NA_circle=utility.fourier_space.coorxy<utility.kmax;

OTF=s_ifft2(s_fft2(utility.NA_circle).*s_fft2(utility.NA_circle));

OTF=OTF./max(OTF(:));
OTF=OTF(floor(size(OTF,1)/2)+(1:size_img(1))-floor(size_img(1)/2),floor(size(OTF,2)/2)+(1:size_img(2))-floor(size_img(2)/2));

sinc_term= sinc(((1:size(OTF,1))'-1-floor(size(OTF,1)/2))./size(OTF,1)).*sinc(((1:size(OTF,2))-1-floor(size(OTF,2)/2))./size(OTF,2));
OTF=OTF.*sinc_term;

end