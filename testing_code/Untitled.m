ii = 170;

figure(100)
subplot(2,3,1)
imagesc(abs(squeeze(Isave(:,:,ii))).^2)
axis image; axis xy
colorbar
subplot(2,3,2)
imagesc(abs(squeeze(Isave(:,:,ii+1))).^2)
axis image; axis xy
colorbar
subplot(2,3,3)
imagesc(abs(squeeze(Isave(:,:,ii+1))).^2-abs(squeeze(Isave(:,:,ii))).^2)
axis image; axis xy
colorbar
subplot(2,3,4)
imagesc(squeeze(ausave(:,:,ii)))
axis image; axis xy
colorbar
subplot(2,3,5)
imagesc(squeeze(ausave(:,:,ii+1)))
axis image; axis xy
colorbar
subplot(2,3,6)
imagesc(squeeze(ausave(:,:,ii+1))-squeeze(ausave(:,:,ii)))
axis image; axis xy
colorbar