function matt = low_pass_filter_2D(matt,cut)

dim_1_indices=(1:size(matt,1))-floor(size(matt,1)/2);
dim_2_indices=(1:size(matt,2))-floor(size(matt,2)/2);
dim_3_indices=(1:size(matt,3))-floor(size(matt,3)/2);

dim_1_indices=dim_1_indices./max(dim_1_indices(:));
dim_2_indices=dim_2_indices./max(dim_2_indices(:));
dim_3_indices=dim_3_indices./max(dim_3_indices(:));

dim_1_indices=reshape(dim_1_indices,[],1,1);
dim_2_indices=reshape(dim_2_indices,1,[],1);
dim_3_indices=reshape(dim_3_indices,1,1,[]);

matt=fftshift(fftn(ifftshift(matt)));

matt=matt.*((dim_1_indices.^2+dim_2_indices.^2+dim_3_indices.^2)<(cut)^2);

matt=fftshift(ifftn(ifftshift(matt)));

end