function OTF=make_bright_field_OTFv2(size_img,resolution,NA,wavelength,cintilator_thickness,cintilator_RI,substrate_thickness,substrate_RI,RI_immersion)


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

utility.k0_cintilator=utility.k0.*cintilator_RI;
utility.k3_cintilator=(utility.k0_cintilator).^2-(utility.fourier_space.coorxy).^2;utility.k3_cintilator(utility.k3_cintilator<0)=0;utility.k3_cintilator=sqrt(utility.k3_cintilator);
utility.refocusing_kernel_cintilator=1i*2*pi*utility.k3_cintilator;
%% spherical aberation

utility.k0_substrate=utility.k0.*substrate_RI;
utility.k3_substrate=(utility.k0_substrate).^2-(utility.fourier_space.coorxy).^2;utility.k3_substrate(utility.k3_substrate<0)=0;utility.k3_substrate=sqrt(utility.k3_substrate);
utility.refocusing_kernel_substrate=1i*2*pi*utility.k3_substrate;

utility.k0_immersion=utility.k0.*RI_immersion;
utility.k3_immersion=(utility.k0_immersion).^2-(utility.fourier_space.coorxy).^2;utility.k3_immersion(utility.k3_immersion<0)=0;utility.k3_immersion=sqrt(utility.k3_immersion);
utility.refocusing_kernel_immersion=1i*2*pi*utility.k3_immersion;

utility.aberrations=exp(utility.refocusing_kernel_substrate*substrate_thickness-utility.refocusing_kernel_immersion*substrate_thickness*RI_immersion/substrate_RI);

%figure; imagesc(utility.NA_circle.*angle(utility.aberrations));axis image;
%%
OTF=0;
for ii=-cintilator_thickness/2:(cintilator_thickness/10):cintilator_thickness/2
OTF=OTF+(abs(s_fft2(utility.NA_circle.*exp(utility.aberrations+utility.refocusing_kernel_cintilator.*ii))).^2);
end
%figure; imagesc(OTF);
OTF=real(s_ifft2(OTF));
%figure; imagesc(real(OTF));
%%ii=10;
%size(utility.NA_circle)
%size(utility.refocusing_kernel)
%figure; imagesc(angle(utility.NA_circle.*exp(utility.refocusing_kernel_cintilator.*ii)))

OTF=OTF./max(OTF(:));
OTF=OTF(floor(size(OTF,1)/2)+(1:size_img(1))-floor(size_img(1)/2),floor(size(OTF,2)/2)+(1:size_img(2))-floor(size_img(2)/2));

sinc_term= sinc(((1:size(OTF,1))'-1-floor(size(OTF,1)/2))./size(OTF,1)).*sinc(((1:size(OTF,2))-1-floor(size(OTF,2)/2))./size(OTF,2));
OTF=OTF.*sinc_term;

end