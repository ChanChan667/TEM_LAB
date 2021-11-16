%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2021  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADF-STEM sample 1: SrTiO3 100
clc;
close all;
clear all;
%% specimen preparation:
lattConst = [3.9051, 3.9051, 0]; % [a b]
sliceDist = [1.9525, 1.9525]; % distance between each slice
expanNum = [5,8]; % expand the unit cell by Expan_Nx = 3 and Expan_Ny = 2, adaptive
% Laters: Each column for an atom
sliceA = [38,   8;...
          1,    1;...
          0,    0.5;...
          0,    0.5];
sliceB = [22,   8,      8;...
          1,    1,      1;...
          0.5,  0.5,    0;...
          0.5,  0,      0.5];

%% sampling:
Lx = expanNum(1) * lattConst(1);
Ly = expanNum(2) * lattConst(2);
Nx = 512;
Ny = 512;
dx = Lx / Nx;
dy = Ly / Ny;
%% STEM settings:
params.KeV = 200;
interCoeff = InteractionCoefficient(params.KeV);
wavLen = HighEnergyWavLen_X(params.KeV);
params.aperture = CircApert_X(Lx, Ly, Nx, Ny, wavLen, 22);
params.Cs3 = 0;
params.Cs5 = 0;
params.df = 0;
params.scanx = linspace(0, 3.9051, 40);
params.scany = linspace(0, 3.9051, 40);

% define multiple STEM detectors:
detectorNum = 2;
params.detector = zeros(Ny, Nx, detectorNum);
params.detector(:, :, 1) = AnnularDetector_X(90, 170, wavLen, Lx, Ly, Nx, Ny);
params.detector(:, :, 2) = AnnularDetector_X(11, 22, wavLen, Lx, Ly, Nx, Ny);
%% Transmission functions:
projPotA = MultiProjPot_conv_X(sliceA, expanNum, lattConst, Lx, Ly, Nx, Ny, 1e-5);
figure;
imagesc(projPotA);
colormap('gray'); axis square;
title('Proj.Pot. A');

projPotB = MultiProjPot_conv_X(sliceB, expanNum, lattConst, Lx, Ly, Nx, Ny, 1e-5);
figure;
imagesc(projPotB);
colormap('gray'); axis square;
title('Proj.Pot. B');
tic;
tfA = exp(1i * interCoeff * projPotA / 1000);
tfB = exp(1i * interCoeff * projPotB / 1000);
tfA = BandwidthLimit(tfA, Lx, Ly, Nx, Ny, 0.67);
tfB = BandwidthLimit(tfB, Lx, Ly, Nx, Ny, 0.67);
transFuncs(:, :, 1) = tfA;
transFuncs(:, :, 2) = tfB;

%% Imaging section:
stackNum = 20;
cbedOption = 0;

theta_0=0.1;
MFP=10;
SEA=MultiSEobjectfunction(theta_0,params,Nx,Ny,Lx,Ly,sliceA, expanNum, lattConst);
SEB=MultiSEobjectfunction(theta_0,params,Nx,Ny,Lx,Ly,sliceB, expanNum, lattConst);
SEobfun(:,:,1)=SEA;
SEobfun(:,:,2)=SEB;

[stemImg,SEimage_foroneslce,SEimage] = STEM_X_SE(Lx, Ly, params, transFuncs,SEobfun,MFP, sliceDist, stackNum,[1,5,10,20,40], 'reduced',...
    cbedOption);
toc;

repStemImg = repmat(stemImg, 6, 6, 1);
stemImg_1 = mat2gray(repStemImg(:, :, 1));
stemImg_2 = mat2gray(repStemImg(:, :, 2));

figure;
imshowpair(stemImg_1, stemImg_2, 'montage');

figure
imagesc(SEimage_foroneslce(:,:,1));
colormap gray;axis square

figure
imagesc(SEimage_foroneslce(:,:,2));
colormao gray; axis square

figure
imagesc(SEimage)
colormap gray; axis square


