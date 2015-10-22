% close all
clear all

 %% Define parameters of files and image set.
exp_.image_dir           = 'C:\Users\soskinne\Documents\Data\BCM Microscope\2015_09_29_PulsingNLS\';
exp_.name_root           = 'Image_';
exp_.images_to_analyze   = [1:3];
exp_.DAPI_channel        = 4;
exp_.Fluo_channel        = 2;
exp_.image_size          = [1024 1024];

display_processing       = true;

if ~exist(exp_.image_dir); mkdir(exp_.image_dir, '\Data Analysis\Segmentation'); end



%% Count the number of slices for each z-stack
fprintf(1,['\nCounting slices'])
exp_.num_z = [];
folder_contents = dir(exp_.image_dir);
for i_im = 1:size(folder_contents,1)   
    file_name = folder_contents(i_im).name;
    % Only look at .tif's
    if ~isempty(regexp(file_name, '.tif', 'once'))
        
        % Requires image filenames to have (# are numbers):
        % name_root #im_num 'z' #zslice 'c' #ch_num file_extension
        nums_in_filename = regexp(file_name, '[0-9]+', 'match');
        im_num = str2num(nums_in_filename{size(nums_in_filename,2)-2});
        zslice = str2num(nums_in_filename{size(nums_in_filename,2)-1});
        
        % loop repeats over all zslices, saving the last
        exp_.num_z(im_num) = zslice;    
        
    end    
    
end

clear folder_contents nums_in_filename file_name im_num zslice 




%% Process each image declared above
for i_im = exp_.images_to_analyze
    
    fprintf(1,['\nLoading and segmenting image ' num2str(i_im) '.'])
    
    
    %% 1) Inialize masks
    DAPI_mask = zeros(exp_.image_size(1),exp_.image_size(2),exp_.num_z(i_im));
    Fluo_mask = zeros(exp_.image_size(1),exp_.image_size(2),exp_.num_z(i_im));


    %% 2) Initial segmentation: edge detection and morph opening
    for i_z = 1:exp_.num_z(i_im);
        
        % Load zslice
        Fluo_temp = double(imread([exp_.image_dir exp_.name_root num2str(i_im,'%03d') 'z' num2str(i_z,'%02g') 'c' num2str(exp_.Fluo_channel) '.tif']));
        DAPI_temp = double(imread([exp_.image_dir exp_.name_root num2str(i_im,'%03d') 'z' num2str(i_z,'%02g') 'c' num2str(exp_.DAPI_channel) '.tif']));
        
        % For DAPI, calc a version of the Sobel filter
        [g_x, g_y] = gradient(DAPI_temp);
        V_DAPI = sqrt(g_x.^2 + g_y.^2);     
        
        % For DAPI, fill holes and open image
        DAPI_mask(:,:,i_z) = imfill(V_DAPI, 'holes');
        DAPI_mask(:,:,i_z) = imopen(DAPI_mask(:,:,i_z), strel('disk',4));
        
        % For Fluorescence, Gaussian filter and open image
        Fluo_mask(:,:,i_z) = imfilter(Fluo_temp, fspecial('gaussian', 5, 2));        
        Fluo_mask(:,:,i_z) = imopen(Fluo_mask(:,:,i_z), strel('disk',4));
       

    end
    
    clear DAPI_temp Fluo_temp V_DAPI norm se_size g_x g_y
   
    
    
    %% 3) Normalize values [Min Max] -> [0.0 1.0]
    DAPI_mask = mat2gray(DAPI_mask);
    Fluo_mask = mat2gray(Fluo_mask);    
    
    
    
    %% 4) Create binary images
    % 4.1) Threshold DAPI images
    % Thresholding method inpired by McHale et al 2011 Development    
    %  -Scans through threshold values 0->1
    %  -For all contiguous regions, inspects area and circularity
    %  -Any region that passes area and circularity requirement
    for i_z = 1:exp_.num_z(i_im) % For each slice
        I_max   = 256;                       % Number of thresholding steps
        im_temp = DAPI_mask(:,:,i_z);        % Make copy of image that will be marked with recognized nuclei
        bw_out  = zeros(exp_.image_size);    % Intialize output image
        low_th  = 500;                       % Area lowerbound
        high_th = 12000;                     % Area upperbound
        circ_th = 0.75;                      % Circularity lowerbound
        for th_I = (1/I_max):(1/I_max):(I_max-1)/I_max
            bw_temp = imfill(im2bw(im_temp,th_I),'holes');      % threshold the image and fill
            bw_temp = imerode(bw_temp,strel('disk',1));         % erode thresholded and filled image
            bw_prop = regionprops(bw_temp,'Area','Perimeter');  % Find contiguous regions in thresholded image
            bw_area = [bw_prop.Area];                           % Find each region's area and circularity
            bw_peri = [bw_prop.Perimeter];
            bw_circ = 4*pi*bw_area./bw_peri.^2;
            ind_true = find(bw_area >= low_th & ...             % Test each region's area and circularity
                            bw_area <= high_th & ...
                            bw_circ >= circ_th );
            bw_true = ismember(bwlabel(bw_temp),ind_true);      % True/False map of regions that meet requirements
            bw_out  = bw_out | bw_true;                         % Update output image with regions that meet requirements
            im_temp(bw_true) = 0;                               % Mark recognized regions (ignore in further thresholds)
        end         
        DAPI_mask(:,:,i_z) = bw_out;
    end
    
    clear I_max bw_area bw_circ bw_out bw_peri bw_prop bw_temp bw_true circ_th high_th
    clear im_temp ind_true low_th th_I
    
    
    % 4.2) Threshold fluorescence images using graythresh (Otsu's method) 
    % Found that occassionally debris will produce crazy fluorescence, messing up thresholding.
    % So, threshold zlices independently and then take max projection 
    for i_z = 1:exp_.num_z(i_im)        
        Fluo_mask(:,:,i_z) = logical(im2bw(Fluo_mask(:,:,i_z), graythresh(Fluo_mask(:,:,i_z))));    
    end
    Fluo_MaxP = max(Fluo_mask, [], 3);
    
    clear Fluo_mask
    
    
    % Display progress of image processing
    if display_processing==true        
        figure
        patch(isosurface(DAPI_mask,.5),'FaceColor','blue','EdgeColor','none');
        view(3); axis equal
        camlight; lighting phong
    end


    %% 5) Morphological processing
    % Open max projection of fluorescence maxP
    Fluo_MaxP = imopen(Fluo_MaxP, strel('disk', 20));
    
    % Join slices in z
    DAPI_mask = imdilate(DAPI_mask, ones(1,1,7));
    % Remove small volumes in DAPI
    DAPI_cc = bwconncomp(DAPI_mask);
    for i_obj = 1:DAPI_cc.NumObjects
        if length(DAPI_cc.PixelIdxList{i_obj}) < 1000
            DAPI_mask(DAPI_cc.PixelIdxList{i_obj}) = 0;
        end
    end
    
    % Round out nuclei
    DAPI_cc = bwconncomp(DAPI_mask);
    DAPI_rp = regionprops(DAPI_cc, 'SubarrayIdx', 'PixelIdxList');
    [ellipse1] = makeellipse(8,4);
    [ellipse2] = makeellipse(4,2);
    for i_obj = 1:DAPI_cc.NumObjects
        idx = [];
        idx = DAPI_rp(i_obj).SubarrayIdx;                   % Get indicies of cell
        cell_mask = DAPI_mask(idx{:});                      % make temporary cell mask
        DAPI_mask(DAPI_rp(i_obj).PixelIdxList) = 0;         % erase cell from mask
        cell_mask = imdilate(cell_mask, ones(1,1,7));       % Join pixels that are separated by 4 pixels in z
        cell_mask = imerode(cell_mask, ellipse1);           % 3d imopen with ellipse shape
        cell_mask = imdilate(cell_mask, ellipse2);          % 3d imopen with ellipse shape
        DAPI_mask(idx{:}) = DAPI_mask(idx{:}) + cell_mask;  % Add processed cell back in
    end
    clear cell_mask DAPI_rp ellipse1 ellipse2
    
    % Remove small volumes in DAPI
    DAPI_cc = bwconncomp(DAPI_mask);
    for i_obj = 1:DAPI_cc.NumObjects
        if length(DAPI_cc.PixelIdxList{i_obj}) < 3000
            DAPI_mask(DAPI_cc.PixelIdxList{i_obj}) = 0;
        end
    end
    
    % Label 3d matrix
    DAPI_cc   = bwconncomp(DAPI_mask);
    DAPI_mask = labelmatrix(DAPI_cc);
    DAPI_MaxP = max(DAPI_mask, [], 3);

    % Adjust fluorescence map to include all pixels in nucleus map
    Fluo_MaxP = Fluo_MaxP ~=0 | DAPI_MaxP ~=0;

    
    % Display progress of image processing
    if display_processing==true        
        figure
        patch(isosurface(DAPI_mask,.5),'FaceColor','blue','EdgeColor','none');
        view(3); axis equal
        camlight; lighting phong
    end

    
    
    
    %% 6) Save file
    save([exp_.image_dir '\Data Analysis\Segmentation\' exp_.name_root 'seg' num2str(i_im, '%03d') '_1.mat'], 'DAPI_mask', 'Fluo_MaxP')
    
    clear L_mask DAPI_cc DAPI_mask Fluo_MaxP DAPI_MaxP


end



