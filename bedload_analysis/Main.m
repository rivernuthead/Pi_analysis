%% Image trasformations and merging
%Author: Marco Redolfi
%Created on: 6-Nov-2015
%ver 2.7 (07-Nov-2016)
%Previous versions:
%	ver 2.2 (01-Jul-2016)
%	ver 2.3 (23-Nov-2016)
%	ver 2.5 (26-Oct-2016)
%   ver 2.6 (02-Nov-2017)

clear all
close all
clc

addpath('Subroutines')


%% Input parameters
Run='W06_Q05r9';        %If empty the last name saved in output folder (from Pi_camera_control) is used
GCP_set='W60';            %Set of GCP used ('W08','W16','W16L' or 'W60')
px_size=1E-3;             %Final pixel size in [m]
distorsion_corr=true;     %If true the radial distortion is corrected
aspect_ratio=3/2;         %Aspect ratio of the original (uncut) images
lens_U='Nikkor_18-55_VR'; %Lens mounted on the upstream camera
lens_D='Nikkor_18-55_ED'; %Lens mounted on the downstream camera
yR_lim=[-0.33 0.33];      %y-limits [m] for cutting the final image (relative to the GCP_centroid), [-0.35 0.35] for W08, [-0.85 0.85] for W16
Nskip=0;                  %Number of images to skip (i.e. because previously computed)
transf_type='projective'; %Type of image projection
show_img=false;           %Showing final image
int_method='nearest';     %Interpolation method used when reprojecting the images
f=27;                     %Focal length [mm], equivalent to a 35 mm film format
Lt=0.40;                  %Length of the interpolation area [m] in the overlapping region
rand_weight=true;         %If true a random choice of the U/D image is used (with weighted probability distribution on a region on Nt points)
max_timelag=15;           %Max timelag allowed between upstream and downstream images (if exceeded a warning message appears)
WB_between=true;          %If true the mean intensity of each RGB band is set to be the same as in the 1st image
WB_within=true;           %If true the color intensity between the two images is equalized
worldfile=false;          %If true a "*.jgw" world file is built, based on the (aligned) system of reference used to specify GCP coordinats
ref_horizontal=false;     %If true the worldfile refers to an horizontal (rather than sloping) plane
plot_GCP=true;            %Plotting GCP location and sequence
rotate_img=true;          %If true images are rotated by 180° (needed when cameras are oriented so that the right part of the image is upstream)

%Finding name of the last Run (if Run=[]);
if isempty(Run)
    list=dir('./Data/Output_cam_ctrl/Output_*');
    for j=1:length(list);
    	date_list(j)=list(j).datenum;
    end
    [date_list,ind]=sort(date_list);
    Run=list(ind(end)).name(8:end-4);
end

%% Reading GCP location
% In a sloping reference system
%   ->x-axis aligned along the longitudinal direction (rails direction)
%   ->y-axis pointing to the right (as in image coordinates)
read_GCP;


%% Reading images

input_dir=['./0_images/',Run];

list_U=dir([input_dir,'/Upstream/*.jpg'])
list_D=dir([input_dir,'/Downstream/*.jpg'])

Nimg=length(list_U);

fname_U=list_U(1).name;
fname_D=list_D(1).name;

U=imread([input_dir,'/Upstream/',fname_U]);
D=imread([input_dir,'/Downstream/',fname_D]);

% Correction parameters:
r=1;
g=1;    %0.95833;
b=1;    %.91011;

U(:,:,1)=uint8(double(U(:,:,1))*r);
U(:,:,2)=uint8(double(U(:,:,2))*g);
U(:,:,3)=uint8(double(U(:,:,3))*b);

%Image info
Uin=imfinfo([input_dir,'/Upstream/',fname_U]);
Din=imfinfo([input_dir,'/Downstream/',fname_D]);
%Shutter speed and aperture
ExpTime_U(1)=Uin.DigitalCamera.ExposureTime;
FNumb_U(1)  =Uin.DigitalCamera.FNumber;
ExpTime_D(1)=Din.DigitalCamera.ExposureTime;
FNumb_D(1)  =Din.DigitalCamera.FNumber;


mkdir(['1_Fused_images/',Run]);

%% White balance and distortion correction

% Barrel distortion
if distorsion_corr
    U = lens_correct(U,lens_U,int_method,aspect_ratio);
    D = lens_correct(D,lens_D,int_method,aspect_ratio);
end

if rotate_img
	%Rotate image by 180°
	U=fliplr(flipud(U));
	D=fliplr(flipud(D));
end

%% From [m] to [px] units
%NB: No need to scale the focal length (has to be in [mm] because relative to the 35 mm frame)

xGCP   =xGCP/px_size;
yGCP   =yGCP/px_size;
zGCP   =zGCP/px_size;
yR_lim =yR_lim/px_size;
Lt     =Lt/px_size;

%% Images transformation based on GCP position

%Upstream image
[tformU,err_project_U,err_geotrasf_U]=pi_rect(U,xGCP(1:6),yGCP(1:6),-zGCP(1:6),f,GCP_set);
title(['Projection error [m]: ',num2str(err_project_U*px_size,'%.4f'),'; Geotransform error [m]: ',num2str(err_geotrasf_U*px_size,'%.4f')])
print(gcf,['./1_Fused_images/Figures/Upstream_',Run,'.jpg'],'-djpeg','-r300')

%Downstream image
[tformD,err_project_D,err_geotrasf_D]=pi_rect(D,xGCP(5:10),yGCP(5:10),-zGCP(5:10),f,GCP_set);
title(['Projection error [m]: ',num2str(err_project_D*px_size,'%.4f'),'; Geotransform error [m]: ',num2str(err_geotrasf_D*px_size,'%.4f')])
print(gcf,['./1_Fused_images/Figures/Downstream_',Run,'.jpg'],'-djpeg','-r300')

RU = imref2d(size(U));
RD = imref2d(size(D));

tic
disp('Fusing images..')

for j=1+Nskip:Nimg

    tic
    fname_U=list_U(j).name;
    fname_D=list_D(j).name;

    datetime_U=datetime(fname_U(1:end-4),'InputFormat','yyyyMMdd-HHmmss');
    datetime_D=datetime(fname_D(1:end-4),'InputFormat','yyyyMMdd-HHmmss');
    time_lag=etime(datevec(datetime_U),datevec(datetime_D));

    if abs(time_lag)>max_timelag
		disp(fname_U)
		disp(fname_D)
        disp(['Warning: elapsed time between upstream and downstream photos >',num2str(max_timelag),'s'])
    end

	U=imread([input_dir,'/Upstream/',fname_U]);
	D=imread([input_dir,'/Downstream/',fname_D]);

	U(:,:,1)=uint8(double(U(:,:,1))*r);
	U(:,:,2)=uint8(double(U(:,:,2))*g);
	U(:,:,3)=uint8(double(U(:,:,3))*b);

	%Saving images at the antecedend time step
	Uin_prec=Uin;
  	Din_prec=Din;

	Uin=imfinfo([input_dir,'/Upstream/',fname_U]);
	Din=imfinfo([input_dir,'/Downstream/',fname_D]);

	%Saving camera properties
	ExpTime_U(j)=Uin.DigitalCamera.ExposureTime;
	FNumb_U(j)  =Uin.DigitalCamera.FNumber;
	ExpTime_D(j)=Din.DigitalCamera.ExposureTime;
	FNumb_D(j)  =Din.DigitalCamera.FNumber;


	%% White balance, registration, fusion, cropping

	Fuse_UD

	if j==1
		%Saving map indicating the proportion of each image used when fusing (1->Upstream; 0->Downstrean; 0-1->Interpolation)
		imwrite(mask,['./1_Fused_images/',Run,'/Mask.tif']);
		%Saving relevant parameters
		save(['./1_Fused_images/',Run,'/Param.mat'],'px_size','distorsion_corr','WB_within','int_method','transf_type','yR_lim')
    end


    %%  Setting mean intensity of the pixel beloning to each camera the same as in the first image

    White_balance


    %% Show final image

    if show_img
		fig=figure
			imshow(C)
		xlabel('x [px]')
		ylabel('y [px]')
    end


    %% Save fused image

    imwrite(C,['./1_Fused_images/',Run,'/Img',num2str(j,'%04d'),'.jpg']);


    %% Save world file

    if worldfile
        %Building refmat
        %   ->we need to convert from [px] to [m]
        if ref_horizontal
			px_size_h=cos(atan(slope)); %px_size once projected on the horizontal plane
            refmat=[[0 -px_size];[px_size_h 0];[px_size_h*(xUL-1) px_size*(-(yUL-1))]];
        else
			refmat=px_size*[[0 -1];[1 0];[xUL-1 -(yUL-1)]];
        end


		worldfilewrite(refmat,['./1_Fused_images/',Run,'/Img',num2str(j,'%04d'),'.jgw']) %Referencing file

	end

    clear C
	disp(['Img',num2str(j,'%04d'),'.jpg'])
    toc

end
