function [res]=page_time_end(A,B)
res=zeros(size(A,1),size(A,2),size(A,3),size(B,4) ,'like',B);
if size(A,4)~=size(B,3)
    error('size missmatch');
end
for ii=1:size(A,3)
    for jj=1:size(B,4)
        %res(:,:,ii,jj)=sum(squeeze(A(:,:,ii,:)).*squeeze(B(:,:,:,jj)),[3]);
        res(:,:,ii,jj)=sum(reshape(A(:,:,ii,:),size(A,1),size(A,2),[]) .* B(:,:,:,jj) ,3);
    end
end
end
