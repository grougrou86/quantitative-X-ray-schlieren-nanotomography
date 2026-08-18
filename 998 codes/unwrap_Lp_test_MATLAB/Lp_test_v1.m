%% [1] Ghiglia, Dennis C., and Louis A. Romero. "Minimum Lp-norm two-dimensional phase unwrapping." JOSA A 13.10 (1996): 1999-2013.



%% problem definition
maxPhase = 0.1*pi;
hmap = single( rgb2gray(imread('peppers.png')) ); % min:zero
% [XX,YY] = meshgrid(1:xx,1:yy);
% hmap = 0.7*XX+0.3*YY;
phase0 = maxPhase*hmap/max(hmap(:));
E = exp(1i*phase0); % wraping.

[yy,xx] = size(hmap);
NAr = 1/3; % too prevent
NAmask = ~mk_ellipse( NAr*xx, NAr*yy, xx, yy );

E2 = ifft2(fft2(E) .* ifftshift(NAmask));
% imagesc(abs(E2)),axis image; colorbar;
% imagesc(angle(E2)),axis image; colorbar;
phi_w = angle(E2);
phi_w = gpuArray(phi_w);
[phi_uw_iter,U_iter,V_iter,residue0]  = unwrap2_Lp(phi_w, 0);

funcTest = @() unwrap2_Lp(phi_w,0);
gputimeit(funcTest,1)

phi_w_CPU = gather(double(phi_w));
funcTest = @() unwrap2(phi_w_CPU);
timeit(funcTest)

%%% show
ax = subplot(231); imagesc(phase0),axis image; colorbar; climFix = ax.CLim; title('truth')
subplot(232); imagesc(phi_w),axis image; colorbar;  title('wrapped phase')
ax = subplot(233); imagesc(phi_uw_iter),axis image; colorbar; ax.CLim = climFix; title('Lp method')
% subplot(234); imagesc(residue0),axis image; colorbar; title(['Residue # = ',num2str(residueN0)])
subplot(235); imagesc(U_iter),axis image; colorbar; title('U_iter');
subplot(236); imagesc(V_iter),axis image; colorbar; title('V0_iter')
% ax = subplot(233); imagesc(phi_uw_iter),axis image; colorbar; title(['p = ',num2str(p)])
ax = subplot(234); imagesc(unwrap2(gather(double(phi_w)))),axis image; colorbar; ax.CLim = climFix; title('path-following method')
