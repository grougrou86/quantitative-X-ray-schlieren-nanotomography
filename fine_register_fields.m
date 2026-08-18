function rawSino_rec_out =fine_register_fields(rawSino_rec)

%%
scann1=-3:0.5:3;
scann2=-0.5:0.1:0.5;
rot1=0;%-0.5:0.1:0.5;%-2:0.5:2;
rot2=0;%-0.5:0.1:0.5;


crp_for_align=5;%5;%3.5;
align_area=rawSino_rec(1+floor(end/crp_for_align):end-floor(end/crp_for_align),1+floor(end/crp_for_align):end-floor(end/crp_for_align),:);
rawSino_rec_out=rawSino_rec;

thetaInd1 = 1:size(align_area,3)/2;
thetaInd2 = size(align_area,3)/2+1:size(align_area,3);

align_area1 = align_area;
align_area2 = gpuArray(single(flip(align_area1(:,:,thetaInd2),2)));
align_area1 = gpuArray(single(align_area1(:,:,thetaInd1)));

align_area2_out=(align_area2) ;
align_area1_out=(align_area1) ;

coo1=(1:size(align_area1,1))-floor(size(align_area1,1)/2)-1;coo1=coo1./size(align_area1,1);
coo1=gpuArray(single(reshape(coo1,[],1)));
coo2=(1:size(align_area1,2))-floor(size(align_area1,2)/2)-1;coo2=coo2./size(align_area1,2);
coo2=gpuArray(single(reshape(coo2,1,[])));

coo1_big=(1:size(rawSino_rec,1))-floor(size(rawSino_rec,1)/2)-1;coo1_big=coo1_big./size(rawSino_rec,1);
coo1_big=gpuArray(single(reshape(coo1_big,[],1)));
coo2_big=(1:size(rawSino_rec,2))-floor(size(rawSino_rec,2)/2)-1;coo2_big=coo2_big./size(rawSino_rec,2);
coo2_big=gpuArray(single(reshape(coo2_big,1,[])));

filt_size=15;%30
filt_size2=10;%30
filt_freq=(1-exp(-(filt_size*coo2).^2)).*(exp(-(filt_size2*coo2).^2));
align_area2 =real(s_ifft2(s_fft2(align_area2).*filt_freq));
align_area1 =real(s_ifft2(s_fft2(align_area1).*filt_freq));
%figure; plot(filt_freq(:));drawnow;
%figure; sliceViewer(gather(align_area2)); error('top');

filter_border=blackman(size(align_area1,1)).*blackman(size(align_area1,2))';

ii_opt=0;
jj_opt=0;
ann_opt=0;
figure;
for aa=1:size(align_area1,3)
    before=align_area1_out(:,:,aa)+align_area2_out(:,:,aa);
    ii_last=0;
    jj_last=0;
    ann_last=0;
    for pp=1:2
        if pp==1
            scann=scann1;
            rot=rot1;
        else
            scann=scann2;
            rot=rot2;
        end
        cost_min=inf;
        for ann=rot

            align_area1_r=imshearotate(align_area1(:,:,aa),ann_last+ann);
            align_area2_r=imshearotate(align_area2(:,:,aa),-(ann_last+ann));
            %align_area1_r=align_area1_r(1+5:end-5,1+5:end-5,:);
            %align_area2_r=align_area2_r(1+5:end-5,1+5:end-5,:);
            for ii=scann
                for jj=scann
                    align_area1_temp=real(s_ifft2(s_fft2(align_area1_r).*exp(-1i.*2.*pi.*((ii_last+ii)*coo1+(jj_last+jj)*coo2))));
                    align_area2_temp=real(s_ifft2(s_fft2(align_area2_r).*exp(+1i.*2.*pi.*((ii_last+ii)*coo1+(jj_last+jj)*coo2))));
                    %cost=sum(filter_border.*abs(align_area1_temp+align_area2_temp),'all');
                    %%cost=sum((abs(align_area1_temp+align_area2_temp)),'all');
                    %cost=sum(filter_border.*sqrt(abs(align_area1_temp+align_area2_temp)),'all');
                    %cost=sum(filter_border.*(abs(align_area1_temp+align_area2_temp)),'all');
                    plus_term=align_area1_temp+align_area2_temp;plus_term(plus_term<0)=0;
                    cost=sum(filter_border.*(abs(plus_term)),'all');
                    
                    %no_neg=align_area1_temp+align_area2_temp;
                    %no_neg(no_neg<0)=0;
                    %cost=sum(filter_border.*abs(no_neg),'all');
                    if cost<cost_min
                        ii_opt=ii;
                        jj_opt=jj;
                        ann_opt=ann;
                        cost_min=cost;
                    end
                end
            end
        end
        ii_last=ii_last+ii_opt;
        jj_last=jj_last+jj_opt;
        ann_last=ann_last+ann_opt;
        %error('stop')
        align_area1_out(:,:,aa)=real(s_ifft2(s_fft2(imshearotate(align_area1_out(:,:,aa),ann_opt)).*exp(-1i.*2.*pi.*(ii_opt*coo1+jj_opt*coo2))));
        align_area2_out(:,:,aa)=real(s_ifft2(s_fft2(imshearotate(align_area2_out(:,:,aa),-ann_opt)).*exp(+1i.*2.*pi.*(ii_opt*coo1+jj_opt*coo2))));
    end

    rawSino_rec_out(:,:,thetaInd1(aa))=real(s_ifft2(s_fft2(imshearotate(rawSino_rec_out(:,:,thetaInd1(aa)),ann_last)).*exp(-1i.*2.*pi.*(ii_last*coo1_big+jj_last*coo2_big))));
    rawSino_rec_out(:,:,thetaInd2(aa))=real(s_ifft2(s_fft2(imshearotate(rawSino_rec_out(:,:,thetaInd2(aa)),ann_last)).*exp(+1i.*2.*pi.*(ii_last*coo1_big-jj_last*coo2_big))));

    %raw_BG_rec_out(:,:,1)=real(s_ifft2(s_fft2(raw_BG_rec_out(:,:,1)).*exp(-1i.*2.*pi.*(ii_last*coo1_big+jj_last*coo2_big))));

    after=align_area1_out(:,:,aa)+align_area2_out(:,:,aa);
    display([num2str(aa) '/' num2str(size(align_area1,3)) ' ---> optimal : ' num2str(ii_last) '/' num2str(jj_last) '/' num2str(ann_last)]);
    %before;
    imagesc(cat(2,before,after)); axis image; drawnow;
end
%error('align tv after field recon')
