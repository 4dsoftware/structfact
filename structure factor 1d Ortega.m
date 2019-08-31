%1-D structure factorclear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% 1d structure factor of 2D spherical colloid system                      %
%                                                                         %
% Programmed by Daniel Du (March 25 2012)                                 %
%                                                                         %
% Both g(r) and g6(r) functions are calculated 10 particle diameters away %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single figure calculation is used below. For batch process see gbatch.m

% Read the tif file of interst---------------------------------------------
%F=imread('C:\Users\Daniel\Reseach\Proj #3 2D tunable interaction\Data\0529\2.0V_50Hz_1.0wt%_0h_2.tif');
F=imread('C:\Users\Daniel\Reseach\Proj #6 Phase Digram\Data\03172012 week\1.4V 1.32wt% Daniel 6\Solid\Data-1.tif');
%F=imread('C:\Users\Optical\Reseach\Proj #6 Phase Digram\Data\03172012 week\1.4V 1.32wt% Daniel 6\sol3.tif');

tic
% Get the posistion of particles by DIP------------------------------------
Fd = double(F(:,:,1));
[sizx,sizy] = size(Fd);
FF = Fd;
b  = bpass(FF,3,40);             % 40 for 2.8 um, 8 for 1 um
bb = mat2gray(b);
pk = pkfnd(bb,0.5,20);           % 20 for 2.8 um, 7 for 1 um
cnt = cntrd(bb,pk,27,0);         % 27 for 2.8 um, 9 for 1 um

% Set Parameters about bin size--------------------------------------------
dr = 0.1;          % Step size in unit of microns
X = 1;             % Range of excluding part
D  = 2.86;          % D is the diameter of a sphere in microns
d  = 10.7;          % Pixel to micron ratio
M  = length(cnt);           % Total number of spheres
N = 0; % start conunting

% Functions ---------------------------------------------------------------
for i = 1:M
    if cnt(i,1)>50+(X+0.8)*D*d & cnt(i,2)>50+(X+0.8)*D*d & cnt(i,1)<sizy-(50+(X+0.8)*D*d) & cnt(i,2)<sizx-(50+(X+0.8)*D*d) %define ROI, exclude the margins of 50pi
        for j = 1:M 
            if cnt(j,1)>50+(X+0.8)*D*d & cnt(j,2)>50+(X+0.8)*D*d & cnt(j,1)<sizy-(50+(X+0.8)*D*d) & cnt(j,2)<sizx-(50+(X+0.8)*D*d)
                N = N + 1;
                Rx(N) = sqrt( (cnt(i,1)-cnt(j,1))^2 + (cnt(i,2)-cnt(j,2))^2 )/d; % Obtaining all possible distances    
            end
        end
    end
end


