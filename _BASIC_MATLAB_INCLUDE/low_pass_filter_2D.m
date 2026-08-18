function matt = low_pass_filter_2D(matt,cut)

dim_1_indices=(1:size(matt,1))-floor(size(matt,1)/2);
dim_2_indices=(1:size(matt,2))-floor(size(matt,2)/2);

dim_1_indices=dim_1_indices./max(dim_1_indices(:));
dim_2_indices=dim_2_indices./max(dim_2_indices(:));

dim_1_indices=reshape(dim_2_indices,length(dim_1_indices),1);

matt=fftshift(fft2(ifftshift(matt)));

matt=matt.*((dim_1_indices.^2+dim_2_indices.^2)<(cut)^2);

matt=fftshift(ifft2(ifftshift(matt)));

end