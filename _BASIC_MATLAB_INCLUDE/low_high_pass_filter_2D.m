function matt = low_high_pass_filter_2D(matt,cut_max,cut_min,smooth)

dim_1_indices=(1:size(matt,1))-floor(size(matt,1)/2);
dim_2_indices=(1:size(matt,2))-floor(size(matt,2)/2);

dim_1_indices=dim_1_indices./max(dim_1_indices(:));
dim_2_indices=dim_2_indices./max(dim_2_indices(:));

dim_1_indices=reshape(dim_2_indices,length(dim_1_indices),1);

matt=fftshift(fft2(ifftshift(matt)));

%matt=matt.*((dim_1_indices.^2+dim_2_indices.^2)<(cut_max)^2).*((dim_1_indices.^2+dim_2_indices.^2)>(cut_min)^2);
filter=...
1./(1+exp(smooth.*(sqrt((dim_1_indices.^2+dim_2_indices.^2))-(cut_max)))).*...
1./(1+exp(-smooth.*(sqrt((dim_1_indices.^2+dim_2_indices.^2))-(cut_min))));
figure; imagesc(filter);
%figure; imagesc(filter);
matt=matt.*filter;
matt=fftshift(ifft2(ifftshift(matt)));

end