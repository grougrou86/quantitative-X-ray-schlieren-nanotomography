RI_reorganized=permute(RI,[2 3 1]);

norm_RI=real(RI_reorganized)./8e-6;
bin_RI=(real(norm_RI)>0.5);
PSF=s_ifft2(...
    s_fft2(norm_RI)...
    .*conj(s_fft2(bin_RI))./((s_fft2(bin_RI)).^2+1000)...
    );

figure; orthosliceViewer(cat(4,norm_RI,bin_RI,norm_RI));
figure; orthosliceViewer(real(PSF));