%% Extracting active area maps from subsequent images
% Result written as differences maps in png format
%ver 2.4 (11-Dec-2016)
%Author: Marco Redolfi
%Created on: 21-May-2015
%Modified on: 22-Jul-2015; 02-Dec-2015; 01-Jul-2016

% Modified on MAY 2023 by Enrico Pandrin.
% Upstream and downstream saturation correction


function Differencing(varargin)

if length(varargin)==1
	Run=varargin{1}
else
	Run=['W06_Q05r9']; %Selected run
end


%% Input parameters

%Analysis parameters
cf=0.7;           %Correction factor that amplifies areas where mask=0 (downstream) 1.0 for runs q07rgm, q15rgm2, 0.7 for run q20rgm2, 0.8 for q10rgm2
band=6;           %1->Red; 2->Green; 3->Blue; 4->RGB average; 5->Hue; 6->Saturation; 7->Value
equalize=0,       %If true images are equalized using the average image Iavg_filt.png (which needs previously calculated using Averaging.m)
Navg_r=5;        %Size of the filter rows (in [px]) (51 or 31 in the old version)
Navg_c=5;		  %Size of the filter columns (in [px])
Nt_avg=2;         %Number of frames for time average. Even numbers only! (6 in the old version)
t_avg_method=1;   %1->Mobile mean; 2->Mobile median
Nskip=0;          %Number of pictures to skip in the analysis
Fampl=1;          %Factor that multiplies the final difference map (used to reduce the relative round off when storing png images)
worldfile=true;   %If true a the "*.jpw" file (if present) of the original images is copied in a "*.pnw" world file


%% Input check

if rem(Nt_avg,2)~=0
    disp('Error! Please provide an even number for Nt_avg')
    stop
end

%% Reading metadata from txt file

data=importdata(['./Data/Data_run/Run_data.txt']);

ind=find(strcmp(data.textdata(2:end),Run));
if isempty(ind)
	disp('Warning: using default Run data!')
	pause(1)
	ind=find(strcmp(data.textdata(2:end),'DEFAULT'));
end

data_line=data.data(ind,:);
Ndata_line=length(data_line);

t_init=data_line(1); %Time before first picture [min]
dt=data_line(2); %Images interval
if Ndata_line>2; y1       =data_line(3); else; y1       =NaN; end ; %Left limit
if Ndata_line>3; y2       =data_line(4); else; y2       =NaN; end ; %Right limit
if Ndata_line>4; x1       =data_line(5); else; x1       =NaN; end ; %Upper limit
if Ndata_line>5; x2       =data_line(6); else; x2       =NaN; end ; %Lower limit
if Ndata_line>6; angle_rot=data_line(7); else; angle_rot=NaN; end ; %Rotation angle

%% Reading mask

ImagesFolder=['./1_Fused_images/',Run];
mask=imread([ImagesFolder,'/Mask.tif']);
mask=double(mask)/255;

%% Reading images

jpegFiles = dir(strcat(ImagesFolder,'/*.jpg'));

%If not found, is the dir can be selected via GUI
if isempty(jpegFiles)
	ImagesFolder=uigetdir('Select the image folder')
	jpegFiles = dir([ImagesFolder,'/*.jpg']);
end

Nphotos=length(jpegFiles);


%Image sorting based on progressive number
for j=1:Nphotos
	fname=jpegFiles(j).name;
	img_numb(j)=str2num(fname(4:end-4));
end
[B,ind] = sort(img_numb); %Index of sorted photos
jpegFiles = jpegFiles(ind); %Sorting images


%Reading average image
if equalize==1
	photo_avg=imread(['./2_Differences/Output/Average_photos/photo_avg_',Run,'.png']);
	imshow(photo_avg)
	Iavg=double(rgb2gray(photo_avg)); %Average garyscale image
	clear photo_avg
	Iavg=Iavg/(mean(mean(Iavg))); %Average image scaled
end

%Reading, rotating and cropping
photo=imread(strcat(ImagesFolder,'/',jpegFiles(Nskip+1).name));
photo=double(photo);
if equalize==1
	for j=1:3
		photo(:,:,j)=photo(:,:,j)./Iavg;
	end
end

%If limits are NaN the image is not cropped
if isnan(x1), x1=1;             , end
if isnan(x2), x2=size(photo,2); , end
if isnan(y1), y1=1;             , end
if isnan(y2), y2=size(photo,1); , end

%Rotate image (only if the angle is not 0 nor NaN)
if abs(angle_rot)>0 %
	photo=imrotate(photo,angle_rot);
end
photo=photo(y1:y2,x1:x2,:);
mask=mask(y1:y2,x1:x2);


%Band extraction (first image)
if band==4
    I=mean(photo);
elseif band<4
    I=photo(:,:,band);
else
	photo_hsv=rgb2hsv(photo);
	I=photo_hsv(:,:,band-4);
	clear photo_hsv
end
I=I.*(1+(1-mask)*(cf-1));


%Filtering matrix
avg_win=ones(Navg_r,Navg_c); % Kernel dimension in rows x columns

frame_count=0; %Frame counter
mkdir(['./2_Differences/Output/',Run]) %Create directory

%Index of the first photo
fname=jpegFiles(1).name;
img_ind_init=str2num(fname(4:end-4));


stack_point=0; %Stack pointer to store images

%% Cycle on different images
for t= 1:Nphotos-Nskip-1

    tic

    disp(['Time step: ',num2str(t)])

	frame_count=frame_count+1; %Frame in the current video

    fname=jpegFiles(Nskip+t+1).name; %Image name
    img_ind=str2num(fname(4:end-4)); %Image number

	time=((img_ind-img_ind_init)*dt+t_init);

	%Reading image of the current cycle
    photo=imread([ImagesFolder,'/',fname]);
	photo=double(photo);
	if equalize==1 %Equalize images on the basis of the average (grayscale) image
		for j=1:3
			photo(:,:,j)=photo(:,:,j)./Iavg;
		end
	end

	%Rotating image (only if the angle is not 0 nor NaN)
	if abs(angle_rot)>0
		photo=imrotate(photo,angle_rot);
	end
    photo=photo(y1:y2,x1:x2,:);

   	%Band selection
    if band==4
        I2=mean(photo);
	elseif band<4
		I2=photo(:,:,band);
	else %hsv bands
		photo_hsv=rgb2hsv(photo);
		I2=photo_hsv(:,:,band-4);
		clear photo_hsv
	end
	
	% Cf coefficient correction:
	% I2=I2.*(1+(1-mask)*(cf-1));




    %Difference between pictures (for the selected band)
    % diff(:,:)=filter2(avg_win,abs(I2-I),'valid')/(Navg_r*Navg_c);

    diff(:,:)=abs(I2-I);

    % delta saturation correction
    % mask is a 1 (upstream) and 0 (downstream) matrix
    disp('PRE mean up and down:');
    diff_up_mean = mean(diff.*abs(mask-1), 'all')
    diff_dw_mean = mean(diff.*mask, 'all')
    delta_mean = diff_dw_mean-diff_up_mean;

    % Perform the correction
    diff(:,:) = diff(:,:) + delta_mean.*(abs(mask-1))*2;
    
    disp('POST mean up and down:');
    diff_up_mean = mean(diff.*abs(mask-1), 'all')
    diff_dw_mean = mean(diff.*mask, 'all')
    



	%Stack differences matrixes to enable time averaging
	stack_point=stack_point+1;
	if stack_point>Nt_avg
		stack_point=1;
	end
	diff_last(:,:,stack_point)=diff;

	if t_avg_method==1 %Mobile mean
		diff=mean(diff_last,3);	
	else %Mobile median
		diff=median(diff_last,3);
    end
	Nlag=size(diff_last,3)/2;

    %Image for the subsequent cycle
    I=I2;

	diff=diff*Fampl;
	if max(max(diff))>255
		disp('Warning. Difference map exceeds 255!')
		pause
    end


    if rem(Nlag,1)==0

       fname_diff=jpegFiles(Nskip+t+1-Nlag).name; %Image name

        imwrite(diff,['./2_Differences/Output/',Run,'/',fname_diff(1:end-4),'.png'],'Compression','none')

        if worldfile %Copying worldfile to the Differences folder (with ".pnw" extension)
            worldfile_path=[ImagesFolder,'/',fname_diff(1:end-4),'.jgw'];
            if exist(worldfile_path) %If the worldfile exists
                copyfile(worldfile_path,['./2_Differences/Output/',Run,'/',fname_diff(1:end-4),'.pnw'])
            end
        end

    end

    toc

end

save(['./2_Differences/Output/',Run,'/Parameters.mat'],'band','equalize','Navg_r','Navg_c','Nt_avg','t_avg_method','Run','Fampl')
