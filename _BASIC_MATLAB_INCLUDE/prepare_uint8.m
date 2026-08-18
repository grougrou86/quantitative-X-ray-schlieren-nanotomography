function [img_out]=prepare_uint8(img,c_range,c_map)
h_sec=img;

h_sec=h_sec-c_range(1);
h_sec=h_sec/(c_range(2)-c_range(1))*size(c_map,1);
h_sec(h_sec<1)=1;
h_sec(h_sec>size(c_map,1))=size(c_map,1);

img_out = ind2rgb(uint8(h_sec-1), c_map);
end