function img=imresize_fft(img,scale)
padd=2;
curve_strength=1/(padd-1);
%curve_strength

%sz_0=size(img);
sz_0=size(img);
sz=sz_0;sz(1)=sz(1).*padd;sz(2)=sz(2).*padd;

img=s_ifft2(padd_crop_to_fit(s_fft2(img),sz));

coo1=reshape(((1:sz(1))-floor(sz(1)/2)-1),[],1);
coo2=reshape(((1:sz(2))-floor(sz(2)/2)-1),1,[]);

coo_r1=gpuArray(single(coo1)).^2.*curve_strength/(2*sz(1));
coo_r2=gpuArray(single(coo2)).^2.*curve_strength/(2*sz(2));

p_i=2*pi*1i;
lens1=exp(p_i.*(coo_r1+coo_r2));
lens2=exp(-(1/scale).*p_i.*(coo_r1+coo_r2));
refocus1=exp((1/curve_strength^2)*(1-scale).*p_i.*(coo_r1+coo_r2));
refocus2=exp((1/curve_strength^2)*-1/(1/scale+1/(1-scale)).*p_i.*(coo_r1+coo_r2));

img=s_fft2(s_ifft2(s_fft2(img.*lens1).*refocus1).*lens2).*refocus2;

img=padd_crop_to_fit(img,sz_0);
img=s_ifft2(img);

