
function scintOTF = scintOTFgen(NAr,deconv_factor,cropSize)

if isscalar(NAr)
    NAr = [NAr,NAr];
end

acSize = 2*round(2*NAr+1);
NAcircle = ~mk_ellipse(NAr(2),NAr(1),acSize(2),acSize(1));

scintPSF = abs(fft2(ifftshift(NAcircle))).^2;
scintOTF = fftshift(ifft2(scintPSF,'symmetric'));
scintOTF = (scintOTF/max(abs(scintOTF(:)))).^deconv_factor;

if nargin < 3
    scintOTF = msize(scintOTF, acSize);                   
else
    scintOTF = msize(scintOTF, cropSize);                   
end

