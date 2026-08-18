function rawSino_rec_out =fine_register_fields(rawSino_rec,fieldReconinputParams)

thetaInd1 = 1:size(rawSino_rec,3)/2;
thetaInd2 = size(rawSino_rec,3)/2+1:size(rawSino_rec,3);

sino2 = flip(rawSino_rec(:,:,thetaInd2),2);
sino1 = rawSino_rec(:,:,thetaInd1);

crp_for_align=5;%5;%3.5;
%sino2=sino2(1+floor(end/crp_for_align):end-floor(end/crp_for_align),1+floor(end/crp_for_align):end-floor(end/crp_for_align),:);
%sino1=sino1(1+floor(end/crp_for_align):end-floor(end/crp_for_align),1+floor(end/crp_for_align):end-floor(end/crp_for_align),:);

coo1=(1:size(sino1,1))-floor(size(sino1,1)/2)-1;coo1=coo1./size(sino1,1);
coo1=gpuArray(single(reshape(coo1,[],1)));
coo2=(1:size(sino1,2))-floor(size(sino1,2)/2)-1;coo2=coo2./size(sino1,2);
coo2=gpuArray(single(reshape(coo2,1,[])));

scann=-0.5:0.1:0.5;

rawSino_rec_out=zeros(size(sino1,1),size(sino1,2),size(sino1,3),length(scann).^2,'single');
ppii=0;
for ii=scann
    for jj=scann
        sino1_temp=real(s_ifft2(s_fft2(sino1).*exp(-1i.*2.*pi.*((ii)*coo1+(jj)*coo2))));
        sino2_temp=real(s_ifft2(s_fft2(sino2).*exp(+1i.*2.*pi.*((ii)*coo1+(jj)*coo2))));
        [Esino, DPCdata] = field_recon_routine_v17(sino1_temp, sino2_temp,1, fieldReconinputParams{:});
        ppii=ppii+1;
        rawSino_rec_out(:,:,:,ppii)=Esino;
    end
end


%rawSino_rec_out=rawSino_rec;
