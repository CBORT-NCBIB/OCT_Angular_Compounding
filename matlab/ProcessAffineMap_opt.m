function out = ProcessAffineMap_opt(Int,procStruct)

% Script to register tomograms to a reference tomogram.
% Inputs:
%   Int : Bscans acquired at different angles.
%   options: A set of options for the algorithm, those must be a struct
%   refInd: index of the central reference B-scan
%   M1StepAngle: step angle between two consequent Bscans
%   h: estimated distance of rotation point from surface in pixels
%   H: distance of rotation from top of tomogram in pixels
%   L: distance of rotation from left edge in pixels
%   isContinuous: Whether the data is acquired with a continuos sweep or not, defualt is false
%   isCoverslip: Whether the data is acquired with a microscope cover slip
%   or not , default is true
% Copyright Bhaskara Rao Chintada (2023), Martin Villiger(2023)


refInd = procStruct.refInd;
M1StepAngle = procStruct.M1angle;


if ~isfield(procStruct,'isContinuous')% compute reference signal independently for even and odd A-lines
    isContinuous = false;
else
    isContinuous = procStruct.isContinuous;
end

if ~isfield(procStruct,'isCoverslip')% compute reference signal independently for even and odd A-lines
    isCoverslip = true;
else
    isCoverslip = procStruct.isCoverslip;
end

imref = Int(:,:,refInd);

nsample = 1.34;
% updated, after calibrating the axial distance:
dx = 4.3734;
dz = 8.596;

% manually generate a matlab affine transform object
ncoverslip = 1.53;
tcoverslip = 50*ncoverslip;% coverslip thickness in pixels
h = procStruct.h;%distance of rotation point from surface in pixels
OPL = procStruct.H;%distance of rotation from top of tomogram in pixels
L = procStruct.L;%distance of rotation from left edge in pixels
alpharef = 0; %ref angle (where we rotate)
  %% try more angles

  %use matlab optimization to find optimal transform
  offset = 1;
  immoving = cat(3,Int(:,:,1:refInd-1),Int(:,:,refInd+1:end));
  StepAngle = M1StepAngle/180*pi/50e3*1024;

  % better estimate of initial solution

  if(isContinuous)
     AlphaLoc = [refInd-1:-1:1 -1:-1:-(refInd-1)]*offset*StepAngle;
  else
     AlphaLoc = [refInd-1:-1:1 -1:-1:-(refInd-1)]*offset*M1StepAngle/180*pi;% current angle
  end


  options = optimset('Display', 'iter','PlotFcns',@optimplotfval);


  if (isCoverslip)
    
    if(isContinuous)
      
      offsetAngle = 20;
      fun = @(x)-findAffineTransform(10*log10(immoving(:,:,refInd-offsetAngle:offsetAngle:refInd+offsetAngle)),10*log10(imref),makeSampleAffineMappingNew(1,nsample,x(1), x(2), x(3), AlphaLoc(refInd-offsetAngle:offsetAngle:refInd+offsetAngle), alpharef, x(4), x(5), tcoverslip,ncoverslip));
      x0 = [h, OPL, L, dz/8.30, dz]+ [10 100 100 0 0].*rand(1);
      [xopt,fval,exitflag,output] = fminsearch(fun,x0);
      
      x0 = xopt;
      fun = @(x)-findAffineTransform(immoving(:,:,1:offsetAngle:end),imref,makeSampleAffineMappingNew(1,nsample,x(1), x(2), x(3), AlphaLoc(1:offsetAngle:end), alpharef, x(4), x(5), tcoverslip, ncoverslip));
      [xopt,fval,exitflag,output] = fminsearch(fun,x0,options);
      
    else
      fun = @(x)-findAffineTransform(10*log10(immoving(:,:,refInd-1:refInd+1)),10*log10(imref),makeSampleAffineMappingNew(1,nsample,x(1), x(2), x(3), AlphaLoc(refInd-1:refInd+1), alpharef, x(4), x(5), tcoverslip,ncoverslip));
      x0 = [h, OPL, L, dz/4.15, dz]+ [10 100 100 0 0].*rand(1);
      [xopt,fval,exitflag,output] = fminsearch(fun,x0);
      
      x0 = xopt;
      fun = @(x)-findAffineTransform(immoving,imref,makeSampleAffineMappingNew(1,nsample,x(1), x(2), x(3), AlphaLoc, alpharef, x(4), x(5), tcoverslip, ncoverslip));
      [xopt,fval,exitflag,output] = fminsearch(fun,x0,options);
%        [xopt,fval,exitflag,output] = fminsearch(fun,x0);
      
    end

  else
     if(isContinuous)
       
      %nsample = 1.5848;
      offsetAngle = 5;
      fun = @(x)-findAffineTransform(10*log10(immoving(:,:,refInd-offsetAngle:offsetAngle:refInd+offsetAngle)),10*log10(imref),makeSampleAffineMappingNew(1,nsample,x(1), x(2), x(3), AlphaLoc(refInd-offsetAngle:offsetAngle:refInd+offsetAngle), alpharef, x(4), x(5)));
      x0 = [h, OPL, L, dz/4.15, dz]+ [10 100 100 0 0].*rand(1);
      [xopt,fval,exitflag,output] = fminsearch(fun,x0);
      
      x0 = xopt;
      fun = @(x)-findAffineTransform(immoving(:,:,1:offsetAngle:end),imref,makeSampleAffineMappingNew(1,nsample,x(1), x(2), x(3), AlphaLoc(1:offsetAngle:end), alpharef, x(4), x(5)));
      [xopt,fval,exitflag,output] = fminsearch(fun,x0,options);
      
     else

      nsample = 1.66;
      fun = @(x)-findAffineTransform(10*log10(immoving(:,:,refInd-1:refInd+1)),10*log10(imref),makeSampleAffineMappingNew(1, nsample, x(1), x(2), x(3), AlphaLoc(refInd-1:refInd+1), alpharef, x(4), x(5)));
      x0 = [h, OPL, L, dz/4.15, dz]+ [10 100 100 0 0].*rand(1);
      [xopt,fval,exitflag,output] = fminsearch(fun,x0);
      
      x0 = xopt;
      fun = @(x)-findAffineTransform(immoving,imref,makeSampleAffineMappingNew(1, nsample, x(1), x(2), x(3), AlphaLoc, alpharef, x(4), x(5)));
      %[xopt,fval,exitflag,output] = fminsearch(fun,x0,options)
      [xopt,fval,exitflag,output] = fminsearch(fun,x0);
     end

  end

  procStruct.dH = xopt(1);
  procStruct.OPL0 = xopt(2);
  procStruct.L = xopt(3);
  procStruct.nsample = nsample;
  procStruct.dx =xopt(4);
  procStruct.dz = xopt(5);
  procStruct.tcoverslip = tcoverslip;
  procStruct.ncoverslip = ncoverslip;
  procStruct.offset = offset;

%%
% go through entire stack

refInt = Int(:,:,refInd);% sum of both polarization channels
clear StackAligned;
% 

for ind = 1:size(Int,3)
    
%     disp(['Mapping B-scan: ' num2str(ind) ' of ' num2str(size(Int,3))]);
    
        if(isContinuous)
            alphaloc = (refInd-ind)*M1StepAngle/180*pi/50e3*1024;
        else
            alphaloc = (refInd-ind)*M1StepAngle/180*pi;
        end
        
        if (isCoverslip)
            tform = makeSampleAffineMappingNew(1,nsample,xopt(1),xopt(2),xopt(3),alphaloc,alpharef,xopt(4), xopt(5), tcoverslip, ncoverslip);
        else
            tform = makeSampleAffineMappingNew(1,nsample,xopt(1),xopt(2),xopt(3),alphaloc,alpharef,xopt(4), xopt(5));
        end
        imtMartin  = imwarp(Int(:,:,ind),tform ,'OutputView',imref2d(size(imref)));
   
    
    if 0
      figure(2);clf;
      subplot(1,3,1);
      imagesc(10*log10(refInt),[60 100]);
      title('ref');
      subplot(1,3,2);
      imagesc(10*log10(Int(:,:,ind)),[60 100]);
      title('original');
      subplot(1,3,3);
      imagesc(10*log10(imtMartin),[60 100]); colormap gray;
      title('Affine mapping');
    end
    
    StackAligned(:,:,ind) = imtMartin;
    tform1a(ind) = tform;
    
end

% make incoherent means
  if(isContinuous)
     totind = 0;
    for ii = 1:8:min(refInd-1,size(StackAligned,3)-refInd)
      totind = totind + 1;
      ImInc(:,:,totind)  = mean(StackAligned (:,:,refInd + [-ii+1:ii-1]),3);
    end
    ImInc(:,:,totind+1) = mean(StackAligned,3);
    
  else
    totind = 0;
    for ii = 1:1:min(refInd-1,size(StackAligned,3)-refInd)
      totind = totind + 1;
      ImInc(:,:,totind)  = mean(StackAligned (:,:,refInd + [-ii+1:ii-1]),3);
    end
  end
%%

procStruct.tform = tform1a;
out.imCoherentPD = refInt;
out.ImInc = ImInc;
out.stack = StackAligned;
out.xopt = procStruct;

end