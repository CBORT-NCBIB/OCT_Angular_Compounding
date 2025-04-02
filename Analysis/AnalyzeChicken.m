%File purpose, create affine 'volumes'
%This program goes to each volume recorded at each M2 position at different
%M1 angles and pulls the each recorded slice at the same position rotated
%by angle M1. It then affine transforms these volume and saves it as an
%'affine volume'


%***User Set Variables**

%Folder where scan volumes are stored in order of aquistion
addpath(genpath('../matlab')); % Functions folder
folderName = '[p.ChickenJune01_03]';
dataRoot = fullfile('../../Output/',folderName,'Tom_z=(400..1800)_x=(32..2048).bin');
outputFolder = fullfile('../../Output/',folderName);


nZBin = length(400:1800); % Samples in Z in .bin file
nXBin = length(32:2048); % Samples in X in .bin file
numAngles = 61; % Total number of angles
stepAngle = 1; % step angle in degrees
% Load tom1Stokes
skipBscans      = 1;%take every nth slice instead of every slice in a volume

fId1 = fopen(fullfile(dataRoot), 'r'); % Open
int = fread(fId1, 'single=>single'); % Read
int = reshape(int, nZBin, nXBin, numAngles); % Reshape
int = int(:,:,1:skipBscans:end);
fclose(fId1); % Close
clear fId1

%***End User Variables**

%General Variables
readOpt = struct;
readOpt.nFrames = numAngles/skipBscans;%How many bscans in the affine volume?

%setup file loading names


%%
testframe = round(readOpt.nFrames/2);

%Display raw images at each M1 position
figure(1);clf;colormap gray;
subplot(1,3,1);

imagesc(10*log10(int(:,:,testframe)),[60 100]);

%**Apply affine transform to collected bscans at different M1 positions**
procStruct.refInd    = ceil(readOpt.nFrames/2);
procStruct.M1angle   = stepAngle*skipBscans;
procStruct.enableMartin = true;
procStruct.isContinuous = false;

centeringImage = 10*log10(mean(int,3));
bwimg = centeringImage>70;
se = strel('disk',5);
bwimg = imclose(bwimg,se);
[r, c] = find(bwimg == 1);
x = mean(c);
y = mean(r);

procStruct.h = 15;%150%distance of rotation point from surface in pixels
procStruct.H = x;%-150;%-350%distance of rotation from top of tomogram in pixels
procStruct.L = y;%500;%512%distance of rotation from left edge in pixels

out = ProcessAffineMap_opt(int,procStruct);
%reference volume
referenceInt = out.imCoherentPD;
%affine volume
affineMartin = out.ImInc(:,:,end);
% affineInt    = out.ImInc1(:,:,end);
%************

loglims = [50 105];

figure(1);
subplot(1,2,1);
imagesc(10*log10(referenceInt),loglims);
title('Reference');
subplot(1,2,2);
imagesc(10*log10(affineMartin),loglims);
title('Martin Affine Mapping Full');
% subplot(1,3,3);
% imagesc(10*log10(affineInt),loglims);
% title('Matlab imregtform Full');
colormap gray;

%%
% 
% Save averaged output as tiff
MovingAvgStack = (10*log10(out.ImInc)-loglims(1))/range(loglims)*65535;
filenameAveraged = sprintf(['Moving_avg_stack_opt_nFrames=%d_skipBscans=%d_dH=%.4f_OPL=%.4f_L=%.4f_nsample=%.4f'...
  '_dx=%.4f_dz=%.4f_tcoverslip=%.4f_ncoverslip=%.4f_offset=%d'],readOpt.nFrames, skipBscans,out.xopt.dH,out.xopt.OPL0,...
  out.xopt.L, out.xopt.nsample, out.xopt.dx ,out.xopt.dz,out.xopt.tcoverslip, out.xopt.ncoverslip, out.xopt.offset);
saveastiff(uint16(MovingAvgStack),fullfile(outputFolder, strcat(filenameAveraged,'.tiff')),struct('overwrite', true));

% Save aligned output as tiff

filenameAligned = sprintf(['Aligned_stack_opt_nFrames=%d_skipBscans=%d_dH=%.4f_OPL=%.4f_L=%.4f_nsample=%.4f'...
  '_dx=%.4f_dz=%.4f_tcoverslip=%.4f_ncoverslip=%.4f_offset=%d'],readOpt.nFrames, skipBscans,out.xopt.dH,out.xopt.OPL0,...
  out.xopt.L, out.xopt.nsample, out.xopt.dx ,out.xopt.dz,out.xopt.tcoverslip, out.xopt.ncoverslip, out.xopt.offset);
AlignedStack = (10*log10(out.stack)-loglims(1))/range(loglims)*65535;
saveastiff(uint16(AlignedStack),fullfile(outputFolder, strcat(filenameAligned,'.tiff')),struct('overwrite', true));

% save settings as .mat file
filenameSettings= sprintf(['Settings_opt_nFrames=%d_skipBscans=%d_dH=%.4f_OPL=%.4f_L=%.4f_nsample=%.4f'...
  '_dx=%.4f_dz=%.4f_tcoverslip=%.4f_ncoverslip=%.4f_offset=%d'],readOpt.nFrames, skipBscans,out.xopt.dH,out.xopt.OPL0,...
  out.xopt.L, out.xopt.nsample, out.xopt.dx ,out.xopt.dz,out.xopt.tcoverslip, out.xopt.ncoverslip, out.xopt.offset);
settings = out.xopt;
save(fullfile(outputFolder,strcat(filenameSettings,'.mat')),'settings') 
%%

% save as a tiff file
filenameRawdata = sprintf('Tom_z=(%d..%d)_x=(%d..%d)_new', 400, 1800, 32,2048);
AlignedStack = (10*log10(int)-loglims(1))/range(loglims)*65535;
saveastiff(uint16(AlignedStack),fullfile(outputFolder, strcat(filenameRawdata,'.tiff')),struct('overwrite', true));

%%
