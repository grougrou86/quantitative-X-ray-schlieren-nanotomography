function display_2D_3D(mat)

mat=squeeze(mat);
color_mode=0;
if size(mat,length(size(mat)))==3
    color_mode=1;
end
if length(size(mat))-color_mode==2
    imagesc(mat);
else
    orthosliceViewer(mat);
end