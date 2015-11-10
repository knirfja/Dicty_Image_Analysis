function [ LcFull, Cell_mask, DAPI_mask ] = DictyCellRec( exp, DAPI_mask )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Cell outlines by thresholding

DIC_temp =  imread([exp.image_dir exp.name_root num2str(exp.image_id, '%03d') 'z' num2str(6, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Disp_channel) '.tif']);
DIC_temp = imadjust(DIC_temp);

% Threshold image
DIC_thresh = im2bw(DIC_temp, graythresh(DIC_temp));
DIC_thresh = imerode(DIC_thresh,strel('disk',5));
DIC_thresh = imdilate(DIC_thresh,strel('disk',2));
DIC_thresh = imclose(DIC_thresh,strel('disk',15));
%     figure(2); imshow(DIC_thresh)

% Create large outline of entire cell mass
DIC_Outline = imclose(DIC_thresh,strel('disk',35));    
%     DIC_Outline = imdilate(DIC_Outline,strel('disk',10));
DIC_Outline=imfill(DIC_Outline,'holes');    
%     figure(4); imshow(DIC_Outline);

% Segment cells by removing the interior edges
DIC_temp= DIC_Outline & ~DIC_thresh;
DIC_temp=imopen(DIC_temp,strel('disk',25));
DIC_temp=imfill(DIC_temp,'holes');
%     figure(5); imshow(temp3);        

% Crop borders
DIC_temp(1:10,:) = 0;
DIC_temp(:,1:10) = 0;
DIC_temp(1014:1024,:) = 0;
DIC_temp(:,1014:1024) = 0;

DAPI_MaxP = max(DAPI_mask, [], 3);
Cell_mask = DAPI_MaxP | DIC_temp;
%     Cell_mask = Cell_mask | RFP_MaxP;

%     figure(7); imshowpair(mat2gray(Disp),DIC_temp);
%     figure(8); imshowpair(mat2gray(Disp),Cell_mask);

%% Cell Segmentation
% For each cell, count the number of nuclei inside.
% If there are none: erase cell
%               > 1: combine nuclei
% Label the cell mask
Temp = bwlabel(Cell_mask);
for i2 = 1:max(Temp(:))

    cell_pix = Temp==i2;
    num_nucl = sort(unique(nonzeros(DAPI_MaxP(cell_pix))));

    if size(num_nucl,1) == 0 
        Cell_mask(cell_pix) = 0;
    elseif size(num_nucl,1) > 1 
        
        
        
        for i3=num_nucl'
            DAPI_mask(DAPI_mask==i3)=min(num_nucl);            
        end
        
%         % OLD segment cells by watershed
%         D = bwdist(nucl_pix);
%         W = int8(watershed(D));
%         W(~cell_pix) = 0;
%         Cell_mask(cell_pix) = -1.*W(cell_pix);
    end
end
% Cell_mask = logical(Cell_mask);

clear cell_pix nucl_pix D W 

%  Match label number of cells to label number of nuclei
LcFull = zeros(size(Cell_mask));
for i2 = int8(unique(nonzeros(DAPI_mask(:))))'
    cell_num = sort(unique(nonzeros(Temp(DAPI_MaxP==i2))));
    LcFull(Temp==cell_num) = -1*i2;
end
LcFull = -1*LcFull;
clear i2 cell_num num_nucl Temp

end

