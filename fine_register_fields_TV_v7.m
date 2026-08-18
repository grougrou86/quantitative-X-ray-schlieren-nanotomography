function rawSino =fine_register_fields_TV_v5(rawSino,figShow)

%rawSino_rec_out_2 =fine_register_fields_TV_v3(rawSino_rec_out,angleLib,fieldReconInputParams);
%[val_numerical grad_numerical]=adiff_numerical(cost, displacement_0, 0.001);[val, grad]=adiff(cost, displacement_0);display(grad);display(grad_numerical);

sino1=padd_crop_to_fit(rawSino,size(rawSino)-[2*round(size(rawSino,1)/6) 2*round(size(rawSino,2)/4) 0]);




Nang = size(sino1,3);

thetaInd1 = 1:Nang/2;
thetaInd2 = Nang/2+1:Nang;

sino2 = flip(sino1(:,:,thetaInd2),2);
sino1 = sino1(:,:,thetaInd1);


displacement_0=zeros(size(sino1,3),2,'single','gpuArray');


sub_cut=3;
for pp=1:sub_cut
    range=round(1+(pp-1).*size(sino1,3)/sub_cut):round(pp.*size(sino1,3)/sub_cut);
    
    temp_sino1=gpuArray(sino1(:,:,range));
    temp_sino2=gpuArray(sino2(:,:,range));
    
    cost=@(displacement) cost_displacement(real(displacement),temp_sino1,temp_sino2);
    %cost(displacement_0)
    
    t_np=1;
    x_n=displacement_0(range,:);
    itt_num=200;
    val_history=zeros(itt_num,1);
    if figShow
        figure;
    end
    for itt=1:itt_num
        itt
        [val, grad]=adiff(cost, displacement_0(range,:));
        val_history(itt)=val;
        [new_displacement_0,x_n,t_np] = FISTA(+grad.*1e9,displacement_0(range,:),x_n,t_np);
        displacement_0(range,:)=new_displacement_0;
        if figShow
            subplot(1,2,1);plot(val_history(1:itt));
            subplot(1,2,2);plot(displacement_0(:,1)); hold on;plot(displacement_0(:,2),'r');hold off;
            drawnow;
        end
    end
end

coo=make_coo(size(rawSino));
coo1=gpuArray(single(coo{1}));
coo2=gpuArray(single(coo{2}));
d1=reshape(displacement_0(:,1),[1 1 size(displacement_0,1)]);
d2=reshape(displacement_0(:,2),[1 1 size(displacement_0,1)]);
trans_filter=reusable(exp(-1i*(2*pi)*(d1.*coo1+d2.*coo2)));
trans_filter2=reusable(exp(-1i*(2*pi)*(-d1.*coo1+d2.*coo2)));
trans_filter=cat(3,trans_filter,trans_filter2);
for kk=1:size(rawSino,3)
    %rawSino(:,:,kk)=sn_fft2(rawSino(:,:,kk));
    %rawSino(:,:,kk)=trans_filter(:,:,kk).*rawSino(:,:,kk);
    %rawSino(:,:,kk)=real(sn_ifft2(rawSino(:,:,kk)));
    
    temp=gpuArray(rawSino(:,:,kk));
    temp=sn_fft2(temp);
    temp=gpuArray(trans_filter(:,:,kk)).*temp;
    rawSino(:,:,kk)=gather(real(sn_ifft2(temp)));
end
end

function cost=cost_displacement(displacement,sino1,sino2)

d1=reshape(displacement(:,1),[1 1 size(displacement,1)]);
d2=reshape(displacement(:,2),[1 1 size(displacement,1)]);


coo_tomo=make_coo(size(sino1));
coo_tomo_1=gpuArray(single(coo_tomo{1}));
coo_tomo_2=gpuArray(single(coo_tomo{2}));
coo_tomo_r2d=sqrt(coo_tomo_1.^2+coo_tomo_2.^2);

%filter=(abs(coo_tomo_r2d)<0.2)&(abs(coo_tomo_r2d)>0.05);
%filter=(abs(coo_tomo_r2d)<0.1)&(abs(coo_tomo_r2d)>0.025);
filter=(abs(coo_tomo_r2d)<0.1).*coo_tomo_r2d.^2;

%filter=(abs(coo_tomo_r2d)<0.2);

sino1=sn_fft2(sino1);sino2=sn_fft2(sino2);
sino1=sino1.*filter;
sino2=sino2.*filter;
sino1=sn_ifft2(sino1);sino2=sn_ifft2(sino2);

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

%cost=mean((abs(sino1-sino2)),'all');
%absorb_val=(real(sino1-sino2));
%cost=mean((abs(sino1-sino2).^2),'all');
cost=mean((abs2(sino1-sino2)),'all');

%cost=mean((abs(sino1+sino2)),'all');
%cost=mean(abs(sino1+sino2)+abs(sino1-sino2),'all');
%cost=mean(abs2(sino1+sino2),'all');

end

