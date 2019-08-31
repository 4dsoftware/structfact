%Status 
%                           Working 
%Comment: No
clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% 2d structure factor and density-density correlation function            %
%             of 2D spherical colloid system                              % 
%                                                                         %
% Programmed by Daniel Du (March 28 2012)                                 %
%                                                                         %
% Both g(r) and g6(r) functions are calculated 10 particle diameters away %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the tif file of interst---------------------------------------------
%F=imread('C:\Users\Daniel\Reseach\Proj #3 2D tunable interaction\Data\0529\2.0V_50Hz_1.0wt%_0h_2.tif');
F=imread('C:\Users\Optical\Reseach\Proj #3 2D tunable interaction\Data\3.8\Gas6\N6.tif');
%F=imread('C:\Users\Optical\Reseach\Proj #6 Phase Digram\Data\03172012 week\1.4V 1.32wt% Daniel 6\sol3.tif');

% Function-----------------------------------------------------------------
Fd = double(F(:,:,1));
[sizx,sizy] = size(Fd);
b  = bpass(Fd,3,40);    %Gaussian filter
sf = real(fftshift(fft2(Fd)));
imagesc (sf);caxis([0 8e5]);
axis([440 960 320 740]);
a = findall(gcf,'type','axes');
set(a,'YTickLabel',[])
set(a,'XTickLabel',[])

% Find where the peaks are-------------------------------------------------
bb0 = mat2gray(sf);
pk0 = pkfnd(bb0,0.15,20);           % 20 for 2.8 um, 7 for 1 um
cnt0 = cntrd(bb0,pk0,27,0);         % 27 for 2.8 um, 9 for 1 um

% Calculate reciprocal lattice vector--------------------------------------
center = cnt0(4,:);
G1 = cnt0(3,:)-cnt0(1,:);
G2 = cnt0(7,:)-cnt0(5,:);
G3 = cnt0(6,:)-center;
G4 = center-cnt0(2,:);
Gvec = (G1+G2+G3+G4)/4;

% Continue to get the DIP data for original image--------------------------
bb = mat2gray(b);
pk = pkfnd(bb,0.5,20);           % 20 for 2.8 um, 7 for 1 um
cnt = cntrd(bb,pk,27,0);         % 27 for 2.8 um, 9 for 1 um
% which includes displacement field information

% Set Parameters about bin size--------------------------------------------
dr = 0.1;          % Step size in unit of microns
X = 10;             % Range of excluding part
D  = 2.86;          % D is the diameter of a sphere in microns
d  = 10.7;          % Pixel to micron ratio
r  = D:dr:X*D;      % Distance variable
M  = length(cnt);           % Total number of spheres
% Function-----------------------------------------------------------------

for nr = 1:length(r)
    %Initialization
    Nc = 0;
    gG_Now = 0; 
    for i = 1:M
        if cnt(i,1)>50+(X+0.8)*D*d & cnt(i,2)>50+(X+0.8)*D*d & cnt(i,1)<sizy-(50+(X+0.8)*D*d) & cnt(i,2)<sizx-(50+(X+0.8)*D*d) %define ROI, exclude the margins of 50pi
            for j =1:M
                R = sqrt( (cnt(i,1)-cnt(j,1))^2 + (cnt(i,2)-cnt(j,2))^2 )/d;
                if R > r(nr) & R <= r(nr)+dr  % Distance within the [r r+dr] range
                    Nc = Nc+ 1;
                    gG_Now = gG_Now + real(conj(exp(1i*(Gvec(1)*cnt(i,1)+Gvec(2)*cnt(i,2))))*exp(1i*(Gvec(1)*cnt(j,1)+Gvec(2)*cnt(j,2))));
                end
            end
        end
    end
    gG(nr) = gG_Now/Nc;
end

plot(r,gG)
