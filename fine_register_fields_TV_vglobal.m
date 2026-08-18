function rawSino =fine_register_fields_TV_vglobal(rawSino,angles,fieldReconinputParams)

fieldReconinputParams_temp=fieldReconinputParams;
fieldReconinputParams_temp{1}=1;
fieldReconinputParams_temp{4}=false;

%rawSino_rec_out_2 =fine_register_fields_TV_v3(rawSino_rec_out,angleLib,fieldReconInputParams);
%[val_numerical grad_numerical]=adiff_numerical(cost, displacement_0, 0.001);[val, grad]=adiff(cost, displacement_0);display(grad);display(grad_numerical);

sino1=padd_crop_to_fit(rawSino,size(rawSino)-[2*round(size(rawSino,1)/6) 2*round(size(rawSino,2)/4) 0]);




Nang = size(sino1,3);

thetaInd1 = 1:Nang/2;
thetaInd2 = Nang/2+1:Nang;

sino2 = flip(sino1(:,:,thetaInd2),2);
sino1 = sino1(:,:,thetaInd1);
angles=angles(thetaInd1);


cost=@(displacement) cost_displacement(real(displacement),sino1,sino2,angles,fieldReconinputParams_temp);
%displacement_0=zeros(size(sino1,3),2,'single','gpuArray');
displacement_0=zeros(1,2,'single','gpuArray');
cost(displacement_0)

t_np=1;
x_n=displacement_0;
itt_num=30;
val_history=zeros(itt_num,1);
figure;
for itt=1:itt_num
    [val, grad]=adiff(cost, displacement_0);
    val_history(itt)=val;
    [displacement_0,x_n,t_np] = FISTA(-grad*10,displacement_0,x_n,t_np);
    subplot(1,2,1);plot(val_history(1:itt));
    subplot(1,2,2);plot(displacement_0(:,1)); hold on;plot(displacement_0(:,2));hold off;
    displacement_0(:,1)%=0
    displacement_0(:,2)%=0.3
    drawnow;
end

coo=make_coo(size(rawSino));
coo1=gpuArray(single(coo{1}));
coo2=gpuArray(single(coo{2}));
d1=reshape(displacement_0(:,1),[1 1 size(displacement_0,1)]);
d2=reshape(displacement_0(:,2),[1 1 size(displacement_0,1)]);
d1=repmat(d1,[1 1 size(sino1,3)]);
d2=repmat(d2,[1 1 size(sino2,3)]);

trans_filter=reusable(exp(-1i*(2*pi)*(d1.*coo1+d2.*coo2)));
trans_filter2=reusable(exp(-1i*(2*pi)*(-d1.*coo1+d2.*coo2)));
trans_filter=cat(3,trans_filter,trans_filter2);

rawSino=sn_fft2(rawSino);
rawSino=trans_filter.*rawSino;
rawSino=real(sn_ifft2(rawSino));
end

function cost=cost_displacement(displacement,sino1,sino2,angles,fieldReconinputParams_temp)

d1=reshape(displacement(:,1),[1 1 size(displacement,1)]);
d2=reshape(displacement(:,2),[1 1 size(displacement,1)]);

coo_tomo=make_coo(size(sino1));
coo_tomo_1=gpuArray(single(coo_tomo{1}));
coo_tomo_2=gpuArray(single(coo_tomo{2}));
coo_tomo_r2d=sqrt(coo_tomo_1.^2+coo_tomo_2.^2);
%{
filter=(abs(coo_tomo_r2d)<0.2)&(abs(coo_tomo_r2d)>0.05);
sino1=sn_fft2(sino1);sino2=sn_fft2(sino2);
sino1=sino1.*filter;
sino2=sino2.*filter;
sino1=sn_ifft2(sino1);sino2=sn_ifft2(sino2);
%}
sino1=padd_crop_to_fit(sino1,size(sino1,[1 2])-[10 10]);
sino2=padd_crop_to_fit(sino2,size(sino2,[1 2])-[10 10]);

coo=make_coo(size(sino1));
coo1=gpuArray(single(coo{1}));coo2=gpuArray(single(coo{2}));

filter2=reusable(exp(-1i*(2*pi)*(d1.*coo1+d2.*coo2)));

sino1=sn_fft2(sino1);sino2=sn_fft2(sino2);
sino1=filter2.*sino1;
sino2=conj(filter2).*sino2;
sino1=sn_ifft2(sino1);sino2=sn_ifft2(sino2);

sino1=reusable(padd_crop_to_fit(sino1,size(sino1,[1 2])-[10 10]));
sino2=reusable(padd_crop_to_fit(sino2,size(sino2,[1 2])-[10 10]));

over_padd=1;

[Esino, DPCdata] = field_recon_routine_v19(real(sino1),real(sino2),1, fieldReconinputParams_temp{:});
 tomo=iradon_adiff(permute(Esino,[2 1 3]),angles,over_padd,true);

%tomo=iradon_adiff(permute((sino1-sino2)+1i.*(sino1+sino2),[2 1 3]),angles,over_padd,true);
%tomo=real(Esino);
sz_crp=floor(size(tomo)./[sqrt(2) sqrt(2) 1]);

tomo=padd_crop_to_fit(tomo,sz_crp);



coo_tomo=make_coo(size(tomo));
coo_tomo_1=gpuArray(single(coo_tomo{1}));
coo_tomo_2=gpuArray(single(coo_tomo{2}));
coo_tomo_3=gpuArray(single(coo_tomo{3}));
tomo=reusable(sn_fftn(tomo));
dx=reusable(sn_ifftn(tomo.*coo_tomo_1));
dy=reusable(sn_ifftn(tomo.*coo_tomo_2));
dz=reusable(sn_ifftn(tomo.*coo_tomo_3));
%TV=sqrt(abs2(real(dx))+abs2(real(dy))+abs2(real(dz)));
TV=sqrt(abs2(imag(dx))+abs2(imag(dy))+abs2(imag(dz)));
%TV=(abs2(imag(dx))+abs2(imag(dy))+abs2(imag(dz)));

cost=mean(TV,'all');
%cost=mean(abs(sino1+sino2)+abs(sino1-sino2),'all');
%cost=mean(abs2(sino1+sino2),'all');

end

