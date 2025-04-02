function [tform,H,h] = makeSampleAffineMappingNew(n1,n2,dH,OPL0,L,Alpha,alpharef,dx,dz,tcoverslip,ncoverslip)
% makes the affine mapping between tomograms recorded in a refractive index
% of n2, behind a planar interface separated from medium n1, at two angular
% positions, rotating the sample including the interface from angle alpha1
% to alpha2, where the rotation point is located H from the 'top' (zero
% pathlength difference), and h from the interface, and L from the left
% side. Rather than directly parameterizing H, and h, we defined dH, the
% distance of the zero-pathlength above the sample surface (not the
% coverslip surface) in addition to OPL0, the pathlength of the apparent
% rotation center.
% dx, dz are the effective sampling distances, in um, of the lateral and
% axial pixel size.

% Includes corrected term for coverslip;

if nargin<10 || isempty(tcoverslip)
    tcoverslip = 0;
end
if nargin<11 || isempty(ncoverslip)
    ncoverslip = 1.53;
end

h = (OPL0-dH*n1-tcoverslip*(ncoverslip-n1))*n1/n2^2;
H = dH + h;

for ind = 1:numel(Alpha)%assuming that alpha may be vectorized
    %T1 and T2 are the affine transforms from sample coordinates (lateral, 
    %axial position) to tomogram coordinates (lateral position, pathlength)
    alpha = Alpha(ind);
    beta = asin(n1/n2*sin(alpha));
    betaslip = asin(n1/ncoverslip*sin(alpha));
    
    T1 = [cos(alpha) -cos(alpha)*tan(beta)*dz/dx h*sin(alpha)*dz/dx + L + tcoverslip*(tan(alpha)-tan(betaslip))*cos(alpha)*dz/dx;...
         n1*sin(alpha)*dx/dz (-n1*tan(beta)*sin(alpha)+n2/cos(beta)) (H*n1-h*n1*cos(alpha)) + ...
         tcoverslip*(ncoverslip/cos(betaslip)-n1/cos(alpha)+(tan(alpha)-tan(betaslip))*sin(alpha)*n1);...
         0,0,1];
    betaref = asin(n1/n2*sin(alpharef));
    T2 = [cos(alpharef) -cos(alpharef)*tan(betaref)*dz/dx h*sin(alpharef)*dz/dx + L + tcoverslip*(tan(alpharef)-tan(betaref))*cos(alpharef)*dz/dx;...
         n1*sin(alpharef)*dx/dz (-n1*tan(betaref)*sin(alpharef)+n2/cos(betaref)) (H*n1-h*n1*cos(alpharef)) + ...
         tcoverslip*(ncoverslip/cos(betaref)-n1/cos(alpharef)+(tan(alpharef)-tan(betaref))*sin(alpharef)*n1);...
         0,0,1];

    TT = T2/T1;
    TT(3,:) = [0,0,1];
    tform(ind) = affine2d(TT');
end