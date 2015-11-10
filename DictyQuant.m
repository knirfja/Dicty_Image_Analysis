% clear all
close all

load_masks=1;
exp.root_dir            = 'C:\Users\frink\Documents\MATLAB\Dicty\2015_09_29_PulsingNLS\'; 
exp.image_dir           = [exp.root_dir 'Images\'];
% exp.image_dir = 'C:\Users\soskinne\Documents\Data\BCM Microscope\2015_09_29_PulsingNLS\';
exp.name_root = 'Image_';
images_to_analyze = 1:107;

exp.dig_z=2;
exp.image_size          = 1024;
exp.image_z     = 15;

exp.DAPI_channel = 4;
exp.GFP_channel = 3;
exp.Fluo_channel = 2;
exp.Disp_channel = 1;

%%
running=true;
i1=1;
g_cellnum=0;
while running
    
    image_id = images_to_analyze(i1);
    fprintf(1,['Quantifying ' num2str(image_id) '\n'])        
    load([exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_2.mat'], 'DAPI_mask', 'LcFull', 'Ann_mask')     
    DAPI_MaxP = max(DAPI_mask, [], 3);
    DAPI_mask=repmat(DAPI_MaxP,[1,1,exp.image_z]);
    LcFull=repmat(LcFull,[1,1,exp.image_z]);
    Ann_mask=repmat(Ann_mask,[1,1,exp.image_z]);
    
    for iter=1:exp.image_z
        DIC(:,:,iter) = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(iter, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Disp_channel) '.tif'])));
        RFP(:,:,iter) = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(iter, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Fluo_channel) '.tif'])));
        DAPI(:,:,iter) = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(iter, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.DAPI_channel) '.tif'])));
        GFP(:,:,iter)  = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(iter, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.GFP_channel) '.tif'])));
    end    
    
    bkgd_GFP=median(GFP(~logical(LcFull)));       
    GFP=GFP-bkgd_GFP;
    GFP(GFP<=0)=0;
    
    bkgd_RFP=median(RFP(~logical(LcFull)));
    RFP=RFP-bkgd_RFP;
    RFP(RFP<=0)=0;
    
    num_cells=unique(nonzeros(DAPI_mask(:,:,:)))';
    for cell_num=num_cells
        % Only include Annuli that are bigger than the nucleus        
        if numel(find(DAPI_mask==cell_num)) <= numel(find(Ann_mask==cell_num))
            g_cellnum=g_cellnum+1;
            for z=1:exp.image_z            
            GFP_temp=GFP(:,:,z);
            
            cell(g_cellnum).nucGFPz(z)=mean(GFP_temp(DAPI_mask(:,:,z)==cell_num));
            cell(g_cellnum).annGFPz(z)=mean(GFP_temp(Ann_mask(:,:,z)==cell_num));
            
%             cell(g_cellnum).nuc_GFP_total=sum(GFP(DAPI_mask==cell_num));
%             cell(g_cellnum).nuc_GFP_px=mean(GFP(DAPI_mask==cell_num));
%             cell(g_cellnum).nuc_RFP_total=sum(RFP(DAPI_mask==cell_num));
%             cell(g_cellnum).nuc_RFP_px=mean(RFP(DAPI_mask==cell_num));
% 
%             cell(g_cellnum).ann_GFP_total=sum(GFP(Ann_mask==cell_num));
%             cell(g_cellnum).ann_GFP_px=mean(GFP(Ann_mask==cell_num));
%             cell(g_cellnum).ann_RFP_total=sum(RFP(Ann_mask==cell_num));
%             cell(g_cellnum).ann_RFP_px=mean(RFP(Ann_mask==cell_num)); 
            
%             cell(g_cellnum).nuc2ann_GFP=cell(g_cellnum).nuc_GFP_px ./ cell(g_cellnum).ann_GFP_px;
            cell(g_cellnum).nuc2ann_GFP=max(cell(g_cellnum).nucGFPz(z) ./ cell(g_cellnum).annGFPz(z));
%             cell(g_cellnum).nuc2ann_RFP=cell(g_cellnum).nuc_RFP_px ./ cell(g_cellnum).ann_RFP_px;
            
            
            end
%             if cell(g_cellnum).nuc2ann_GFP==0
%                 disp('')
%             end
%             if isnan(cell(g_cellnum).nuc2ann_GFP)
%                 disp('')
%             end
        end
    end
    clear cell_num        
    
    i1=i1+1;
    if i1 > numel(images_to_analyze)
        running=false;
    end                
end
%%
all_nuc2ann_GFP=[cell.nuc2ann_GFP];
all_nuc2ann_GFP(all_nuc2ann_GFP>4)=4;

% group1_nuc2ann_GFP

figure
h1=histogram(all_nuc2ann_GFP,20);

%% Unused
% %%
% running=true;
% i1=1;
% g_cellnum=0;
% while running
%     
%     image_id = images_to_analyze(i1);
%     fprintf(1,['Quantifying ' num2str(image_id) '\n'])        
%     load([exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_2.mat'], 'DAPI_mask', 'LcFull', 'Ann_mask')     
%     DAPI_MaxP = max(DAPI_mask, [], 3);
%     DAPI_mask=repmat(DAPI_MaxP,[1,1,exp.image_z]);
%     LcFull=repmat(LcFull,[1,1,exp.image_z]);
%     Ann_mask=repmat(Ann_mask,[1,1,exp.image_z]);
%     
%     for iter=1:exp.image_z
%         DIC(:,:,iter) = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(iter, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Disp_channel) '.tif'])));
%         RFP(:,:,iter) = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(iter, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Fluo_channel) '.tif'])));
%         DAPI(:,:,iter) = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(iter, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.DAPI_channel) '.tif'])));
%         GFP(:,:,iter)  = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(iter, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.GFP_channel) '.tif'])));
%     end    
%     
%     bkgd_GFP=median(GFP(~logical(LcFull)));       
%     GFP=GFP-bkgd_GFP;
%     GFP(GFP<=0)=0;
%     
%     bkgd_RFP=median(RFP(~logical(LcFull)));
%     RFP=RFP-bkgd_RFP;
%     RFP(RFP<=0)=0;
%     
%     num_cells=unique(nonzeros(DAPI_mask(:,:,:)))';
%     for cell_num=num_cells
%         % Only include Annuli that are bigger than the nucleus        
%         if numel(find(DAPI_mask==cell_num)) <= numel(find(Ann_mask==cell_num))
%             g_cellnum=g_cellnum+1;
%             
%             
%             
%             cell(g_cellnum).nuc_GFP_total=sum(GFP(DAPI_mask==cell_num));
%             cell(g_cellnum).nuc_GFP_px=mean(GFP(DAPI_mask==cell_num));
%             cell(g_cellnum).nuc_RFP_total=sum(RFP(DAPI_mask==cell_num));
%             cell(g_cellnum).nuc_RFP_px=mean(RFP(DAPI_mask==cell_num));
% 
%             cell(g_cellnum).ann_GFP_total=sum(GFP(Ann_mask==cell_num));
%             cell(g_cellnum).ann_GFP_px=mean(GFP(Ann_mask==cell_num));
%             cell(g_cellnum).ann_RFP_total=sum(RFP(Ann_mask==cell_num));
%             cell(g_cellnum).ann_RFP_px=mean(RFP(Ann_mask==cell_num)); 
%             
%             cell(g_cellnum).nuc2ann_GFP=cell(g_cellnum).nuc_GFP_px ./ cell(g_cellnum).ann_GFP_px;
%             cell(g_cellnum).nuc2ann_RFP=cell(g_cellnum).nuc_RFP_px ./ cell(g_cellnum).ann_RFP_px;
% %             if cell(g_cellnum).nuc2ann_GFP==0
% %                 disp('')
% %             end
% %             if isnan(cell(g_cellnum).nuc2ann_GFP)
% %                 disp('')
% %             end
%         end
%     end
%     clear cell_num        
%     
%     i1=i1+1;
%     if i1 > numel(images_to_analyze)
%         running=false;
%     end                
% end