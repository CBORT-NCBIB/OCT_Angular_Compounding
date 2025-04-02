function [ccsum,cc] = findAffineTransform(Immoving,imref,Tform)

for ind = 1:numel(Tform)
    tform = Tform(ind);
    immoving = Immoving(:,:,ind);
    imcorr = imwarp(immoving,tform,'OutputView',imref2d(size(imref)));
%     mask = (imref>0);
%     cc(ind) = sum(sum(imcorr.*imref.*mask))./sqrt(sum(imcorr(:).^2.*mask(:)))./sqrt(sum(imref(:).^2.*mask(:)));
    cc(ind) = sum(sum(imcorr.*imref))./sqrt(sum(imcorr(:).^2))./sqrt(sum(imref(:).^2));
end
ccsum = sum(cc);
