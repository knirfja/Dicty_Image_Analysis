function Cell_mask = DictyCellRecOLD(exp_,i_im)

fprintf(1,['Loading image and segmenting cell mask ' num2str(i_im) '.\n'])


%% 1) Inialize masks
Phase_mask = zeros(exp_.image_size(1),exp_.image_size(2),exp_.num_z(i_im));
RFP_mask = zeros(exp_.image_size(1),exp_.image_size(2),exp_.num_z(i_im));
GFP_mask = zeros(exp_.image_size(1),exp_.image_size(2),exp_.num_z(i_im));
%     DAPI_org = zeros(exp_.image_size(1),exp_.image_size(2),exp_.num_z(i_im));


%% 2) Initial segmentation: edge detection and morph opening
%     for i_z = 1:exp_.num_z(i_im);
for i_z = 1

    % Load zslice
    RFP_temp = double(imread([exp_.image_dir exp_.name_root num2str(i_im,'%03d') 'z' num2str(i_z,'%02g') 'c' num2str(exp_.RFP_channel) '.tif']));
    GFP_temp = double(imread([exp_.image_dir exp_.name_root num2str(i_im,'%03d') 'z' num2str(i_z,'%02g') 'c' num2str(exp_.GFP_channel) '.tif']));
    Phase_temp = double(imread([exp_.image_dir exp_.name_root num2str(i_im,'%03d') 'z' num2str(i_z,'%02g') 'c' num2str(exp_.Phase_channel) '.tif']));

    % Debug#
    figure; imshow(GFP_temp,[]); figure; imshow(RFP_temp,[]); figure; imshow(Phase_temp,[]);

    % Normalize values [Min Max] -> [0.0 1.0]

    RFP_temp = mat2gray(RFP_temp);    
    GFP_temp = mat2gray(GFP_temp);
    Phase_temp = mat2gray(Phase_temp);       

    % For Fluorescence, Gaussian filter and open image
    RFP_mask(:,:,i_z) = imfilter(RFP_temp, fspecial('gaussian', 5, 2));        
    RFP_mask(:,:,i_z) = imopen(RFP_mask(:,:,i_z), strel('disk',4));

    [pg_x{i_z}, pg_y{i_z}] = gradient(Phase_temp);
    V_Phase = sqrt(pg_x{i_z}.^2 + pg_y{i_z}.^2);   

    [gg_x{i_z}, gg_y{i_z}] = gradient(GFP_temp);
    V_GFP = sqrt(gg_x{i_z}.^2 + gg_y{i_z}.^2);                  

    Phase_mask(:,:,i_z) = imfill(V_Phase, 'holes');
    Phase_mask(:,:,i_z) = imopen(Phase_mask(:,:,i_z), strel('disk',4));   

    % For Fluorescence, Gaussian filter and open image
    RFP_temp = imfilter(RFP_temp, fspecial('gaussian', 5, 2));        
    RFP_temp = imopen(RFP_temp, strel('disk',4));

    [rg_x{i_z}, rg_y{i_z}] = gradient(RFP_temp);
    V_RFP = sqrt(rg_x{i_z}.^2 + rg_y{i_z}.^2);  

    V_PG_composite = sqrt(pg_x{i_z}.^2 + pg_y{i_z}.^2 + gg_x{i_z}.^2 + gg_y{i_z}.^2);
    V_PG_composite2 = sqrt(pg_x{i_z}.^2 + pg_y{i_z}.^2 + gg_x{i_z}.^2 + gg_y{i_z}.^2 + rg_x{i_z}.^2 + rg_y{i_z}.^2);

    RFP_mask(:,:,i_z) = RFP_temp;
end

%     clear RFP_temp GFP_temp Phase_temp


%% 3) Normalize values [Min Max] -> [0.0 1.0]

RFP_mask = mat2gray(RFP_mask);    
GFP_mask = mat2gray(GFP_mask);
Phase_mask = mat2gray(Phase_mask);       


% 4.2) Threshold fluorescence images using graythresh (Otsu's method) 
% Found that occassionally debris will produce crazy fluorescence, messing up thresholding.
% So, threshold zlices independently and then take max projection 
for i_z = 1:exp_.num_z(i_im)        
    RFP_mask(:,:,i_z) = logical(im2bw(RFP_mask(:,:,i_z), graythresh(RFP_mask(:,:,i_z))));    
end
RFP_MaxP = max(RFP_mask, [], 3);        



%% 5) Morphological processing
% Open max projection of fluorescence maxP
RFP_MaxP = imopen(RFP_MaxP, strel('disk', 20));


Cell_mask(i_im) = RFP_MaxP

end