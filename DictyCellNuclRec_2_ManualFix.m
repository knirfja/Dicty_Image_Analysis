% clear all
close all

load_masks=1;
autoprocess=0;
exp.root_dir            = 'C:\Users\frink\Documents\MATLAB\Dicty\2015_09_29_PulsingNLS\'; 
exp.image_dir           = [exp.root_dir 'Images\'];
% exp.image_dir = 'C:\Users\soskinne\Documents\Data\BCM Microscope\2015_09_29_PulsingNLS\';
exp.name_root = 'Image_';
images_to_analyze = 25:107;

exp.dig_z=2;
exp.image_size          = 1024;

exp.DAPI_channel = 4;
exp.GFP_channel = 3;
exp.Fluo_channel = 2;
exp.Disp_channel = 1;
      
%%
running=true;
i1=1;
while running  
    image_id = images_to_analyze(i1);
    exp.image_id=image_id;
        %% Prepare all images
        fprintf(1,['Manually Adjusting Mask ' num2str(image_id) '\n'])
        %load([exp.image_dir '\Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_test_2.mat'], 'DAPI_mask')
        load([exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_1.mat'])
    %     Cell_mask=RFP_MaxP;


        DIC = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(6, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Disp_channel) '.tif'])));
        Fluo = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(6, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Fluo_channel) '.tif'])));
        DAPI = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(6, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.DAPI_channel) '.tif'])));
        GFP  = double(imadjust(imread([exp.image_dir exp.name_root num2str(image_id, '%03d') 'z' num2str(6, ['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.GFP_channel) '.tif'])));

        

        slices = size(DAPI_mask,3);

    if load_masks
        load([exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_2.mat'], 'DAPI_mask', 'LcFull', 'Ann_mask')        
        cc = bwconncomp(logical(DAPI_mask));
        DAPI_MaxP = max(DAPI_mask, [], 3);
        [ Ann_mask ] = CreateAnnulus(DAPI_MaxP, LcFull);  
    else
        %% Cell Rec
        cc = bwconncomp(logical(DAPI_mask));
        DAPI_mask = labelmatrix(cc);
        [ LcFull, Cell_mask, DAPI_mask ] = DictyCellRec( exp, DAPI_mask );
        DAPI_MaxP = max(DAPI_mask, [], 3);
        [ DAPI_mask, DAPI_MaxP] = LabelMatchNuc2Cell( DAPI_mask,DAPI_MaxP,LcFull );
        %% Create Annulus
        [ Ann_mask ] = CreateAnnulus(DAPI_MaxP, LcFull);  
    end
    
    %% Editing        
    editing = true;
    saved_flag=true;
    while editing == true;                
        
        %% Display
        % Make DAPI and DIC composites
        R_DAPI = zeros(size(DAPI));
        G_DAPI = zeros(size(DAPI));
        B_DAPI = DAPI;
        R_DIC = DIC;
        G_DIC = DIC;
        B_DIC = DIC;
        R_GFP = zeros(size(GFP));
        G_GFP = GFP;
        B_GFP = zeros(size(GFP));
                
        Fluo_perim = bwperim(LcFull);
        DAPI_perim = bwperim(DAPI_MaxP);
        R_DIC(Fluo_perim==1)  = 65535;
        R_DAPI(Fluo_perim==1) = 65535;
        R_DAPI(DAPI_perim==1) = 65535;
        G_DAPI(DAPI_perim==1) = 65535;
        B_DAPI(DAPI_perim==1) = 65535;
        
        RGB_DIC_outline  = mat2gray(cat(3,R_DIC, G_DIC, B_DIC));
        RGB_DAPI_outline = mat2gray(cat(3,R_DAPI,G_DAPI,B_DAPI));
        
        clear R G B Fluo_perim DAPI_perim R_DAPI G_DAPI B_DAPI R_DIC G_DIC B_DIC
        
        DIC_Ann_rgb=cat(3,DIC, DIC, DIC);
        Ann_rgb = label2rgb(double(Ann_mask), jet(double(max(Ann_mask(:)))));
        Ann_rgb = im2uint16(Ann_rgb);
        
        
        R_GFP(bwperim(Ann_mask)) = 65535; %Ann_rgb(bwperim(Ann_mask));
        G_GFP(bwperim(Ann_mask)) = 65535; %Ann_rgb(bwperim(Ann_mask));
        B_GFP(bwperim(Ann_mask)) = 65535; %Ann_rgb(bwperim(Ann_mask));                
        RGB_GFP_outline=mat2gray(cat(3,R_GFP,G_GFP,B_GFP));
        
        for iter=1:3
            color_Ann = Ann_rgb(:,:,iter);
            color_DIC2 = DIC_Ann_rgb(:,:,iter);
            
            color_DIC2(logical(bwperim(Ann_mask))) = color_Ann(logical(bwperim(Ann_mask)));
            DIC_Ann_rgb(:,:,iter)=color_DIC2;            
        end    
        
        clear iter color_Ann color_DIC2 Ann_rgb                 

        fig_size=450;
        figure(1)
        imshow(mat2gray(RGB_DAPI_outline))
        set(gcf, 'Position', [5, 600, fig_size, fig_size])

        figure(2)
        imshow(mat2gray(RGB_DIC_outline))
        set(gcf, 'Position', [fig_size+20, 600, fig_size, fig_size])
    
        figure(3)
        imshow(mat2gray(DAPI_MaxP))
        colormap jet
        set(gcf, 'Position', [2*fig_size+35, 600, fig_size, fig_size])
        
        figure(5)
%         imshow(mat2gray(DIC_Ann_rgb))
        imshow(mat2gray(RGB_GFP_outline))
        set(gcf, 'Position', [3*fig_size+50, 600, fig_size, fig_size])                                           

        figure(10)
        imshow(mat2gray(Fluo))
        set(gcf, 'Position', [565, 200, 300, 300])
        
        figure(11)
        imshow(mat2gray(DIC))
        set(gcf, 'Position', [205, 200, 300, 300])
        
        figure(12)
        imshow(mat2gray(GFP))
        set(gcf, 'Position', [925, 200, 300, 300])
        
%         figure(4)
%         patch(isosurface(DAPI_mask,.5),'FaceColor','blue','EdgeColor','none');
%         view(3); axis equal
%         camlight; lighting phong
    
        figure(3)  
        
        if autoprocess
            button = '.';
        else
            %% Input
            [x,y,button] = ginput(1);
        end        
        
        if button == 'q' %% Quit - Done with this frame
            if saved_flag
                editing = false;
                running = false;
            else
                save_button=questdlg('Save before quitting image?');
                if ~strcmp(save_button,'Cancel')
                    editing = false;
                    running = false;
                end  
            end               
        elseif button == '>' || button == '.' %% Go to next image
            if saved_flag
                editing = false;                
                i1=i1+1;
                if i1 > numel(images_to_analyze)
                    running = false;
                end
            else
                save_button=questdlg('Save before quitting image?');
                if strcmp(save_button,'Cancel')
                else
                    editing = false;                
                    i1=i1+1;
                    if i1 > numel(images_to_analyze)
                        running = false;
                    end
                end
            end            
        elseif button == '<' || button == ',' %% Go to next image
            if saved_flag
                editing = false;                
                i1=i1-1;
                if i1 < 1
                    i1=1;
                end
            else
                save_button=questdlg('Save before quitting?');
                if strcmp(save_button,'Cancel')
                else
                    editing = false;                
                    i1=i1-1;
                    if i1 < 1
                        i1=1;
                    end
                end
            end            
            
        elseif button == 'a' %% Add - Perform segmentation on the region around this cell.
            
            figure(4)
            BW = roipoly(mat2gray(DAPI));
            close(4)
            
            slices = size(DAPI_mask,3);
            Fluo_add = zeros(exp.image_size,exp.image_size,slices(end));
            DAPI_add = zeros(exp.image_size,exp.image_size,slices(end));
            
            %%%  Load and detect edges
            for i2 = 1:slices;
                
                Fluo_temp = double(imread([exp.image_dir exp.name_root num2str(image_id,'%03d') 'z' num2str(i2,['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Fluo_channel) '.tif']));
                DAPI_temp = double(imread([exp.image_dir exp.name_root num2str(image_id,'%03d') 'z' num2str(i2,['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.DAPI_channel) '.tif']));
        
                [g_x, g_y] = gradient(DAPI_temp);
                V_DAPI = sqrt(g_x.^2 + g_y.^2);

                Fluo_temp(BW==0) = min(Fluo_temp(:));
                V_DAPI(BW==0) = min(V_DAPI(:));

                Fluo_add(:,:,i2) = imfilter(Fluo_temp, fspecial('gaussian', 5, 2));
                DAPI_add(:,:,i2) = imfill(V_DAPI, 'holes');

                Fluo_add(:,:,i2) = imopen(Fluo_add(:,:,i2), strel('disk',4));
                DAPI_add(:,:,i2) = imopen(DAPI_add(:,:,i2), strel('disk',4));
                

            end

            clear DAPI_temp Fluo_temp V_DAPI norm se_size g_x g_y

            DAPI_add = mat2gray(DAPI_add);
            Fluo_add = mat2gray(Fluo_add);

            threshold_Fluo = 1;

            
            %% Create binary mask
            for i2 = 1:slices   

                im = DAPI_add(:,:,i2);

                %% Heng's method
                I_max = 256;
                im_temp = im;
                bw_im = zeros(size(im));
                low_th = 1000;
                high_th = 12000;
                cir_th = 0.75;
                conv_th = 1.2;
                tic
                for I = 1:(I_max-1)
                    th_I = I/I_max;
                    bw_temp = imfill(im2bw(im_temp,th_I),'holes');
                    bw_temp = imerode(bw_temp,strel('disk',1));
                    bw_prop = regionprops(bw_temp,'Area','Perimeter','ConvexArea');
                    bw_area = [bw_prop.Area];
                    bw_perim = [bw_prop.Perimeter];
                    bw_conv = [bw_prop.ConvexArea];
                    ind_true = find(bw_area >= low_th & bw_area <= high_th & 4*pi*bw_area./bw_perim.^2 >= cir_th & bw_conv./bw_area <= conv_th);
                    bw_true = ismember(bwlabel(bw_temp),ind_true);
                    bw_im = bw_im | bw_true;
                    im_temp(bw_true) = 0; 
                end  
                toc

                clear bw_true im_temp bw_temp im
        %         figure(1)
        %         imagesc(bw_im)

                la_im = bwlabel(bw_im);
                bw_im = zeros(size(bw_im));
                Bkgd_cell = zeros(size(bw_im));
                

                for i3 = 1:max(la_im(:))

                    DAPI_temp = mat2gray((imread([exp.image_dir exp.name_root num2str(image_id,'%03d') 'z' num2str(i2,'%02g') 'c' num2str(exp.DAPI_channel) '.tif'])));

                    Bkgd_pix = DAPI_temp(la_im==0);
                    Bkgd_cell(la_im==0 | la_im==i3) = DAPI_temp(la_im==0 | la_im==i3);
                    Bkgd_cell_pix = DAPI_temp(la_im==0 | la_im==i3);

                    Bkgd_cell_pix_sort = sort(Bkgd_cell_pix);
                    th = Bkgd_cell_pix_sort(round(.9*length(Bkgd_cell_pix_sort)));

                    th_Bkgd_cell = graythresh(Bkgd_cell_pix);

                    bw_DAPI = im2bw(DAPI_temp, th_Bkgd_cell);

                    bw_cell = bw_DAPI & la_im==i3;

        %             figure(1)
        %             imagesc(bw_im)

                    bw_im = bw_im | bw_cell;


        %             im_temp2 = im_temp;
        %             im_temp2(la_im==i3) = im(la_im==i3);
        %             th = graythresh(im_temp2(im_temp2~=0));
        %             cell_bw = im2bw(im_temp2,th);
        %             bw_im = bw_im | ( cell_bw & la_im==i3 );
                end


                DAPI_add(:,:,i2) = bw_im;
                %% Normal thresholding
                %DAPI_mask(:,:,i2) = logical(im2bw(DAPI_mask(:,:,i2), threshold_DAPI));
                Fluo_add(:,:,i2) = logical(im2bw(Fluo_add(:,:,i2), threshold_Fluo));

        %         figure(2)
        %         subplot(1,2,1)
        %         imagesc(DAPI_mask(:,:,i2))
        %         subplot(1,2,2)
        %         imagesc(DAPI_mask_Heng(:,:,i2))

        %       %figure;imagesc(DAPI_mask(:,:,i2))
                %figure;imagesc(Fluo_mask(:,:,i2))
            end

            clear Bkgd_cell Bkgd_cell_pix Bkgd_cell_pix_sort Bkgd_pix bw_DAPI bw_cell
            
            %%  Fix up until here.
            
            new_cell = max(Fluo_add,[],3); %Max projection
            new_nucl = DAPI_add;

            %figure;imagesc(Cell_mask)

            clear Fluo_mask BW_Fluo BW_DAPI DAPI_add Fluo_add threshold_DAPI threshold_Fluo DAPI_temp V norm se_size


            fprintf(1,['\n            Processing image ' num2str(image_id) '.'])
            %%% Morphological processing
            % Open Fluor maxP
            new_cell = imopen(new_cell, strel('disk', 20));

            % Remove small volumes in DAPI
            nucl_cc = bwconncomp(new_nucl);
            for i2 = 1:nucl_cc.NumObjects
                if length(nucl_cc.PixelIdxList{i2}) < 1000
                    new_nucl(nucl_cc.PixelIdxList{i2}) = 0;
                end
            end

            % Round out nuclei
            nucl_cc = bwconncomp(new_nucl);
            nucl_rp = regionprops(nucl_cc, 'SubarrayIdx', 'PixelIdxList');
            [ellipse] = makeellipse(8,4);
            for i2 = 1:nucl_cc.NumObjects
                mask_zfix = [];
                mask_open = [];
                idx = nucl_rp(i2).SubarrayIdx;                      % Get indicies of cell
                cell_mask = new_nucl(idx{:});                      % make temporary cell mask
                new_nucl(nucl_rp(i2).PixelIdxList) = 0;            % erase cell from mask
                cell_mask = imdilate(cell_mask, ones(1,1,5));       % Join pixels that are separated by 4 pixels in z
                cell_mask = imopen(cell_mask, ellipse);             % 3d imopen with ellipse shape
                new_nucl(idx{:}) = new_nucl(idx{:}) + cell_mask;  % Add processed cell back in
            end
            clear cell_mask DAPI_rp ellipse

            % Remove small volumes in DAPI
            nucl_cc = bwconncomp(new_nucl);
            for i2 = 1:nucl_cc.NumObjects
                if length(nucl_cc.PixelIdxList{i2}) < 1000
                    new_nucl(nucl_cc.PixelIdxList{i2}) = 0;
                end
            end

            % Label 3d matrix
            nucl_cc   = bwconncomp(new_nucl);
            new_nucl = labelmatrix(nucl_cc);
            nucl_MaxP = max(new_nucl, [], 3);

            % Fix Fluo maxp
            new_cell = new_cell ~=0 | nucl_MaxP ~=0;
              
            Cell_mask(new_cell==1) = 1;
            DAPI_mask(new_nucl==1) = 1;
            
            cc = bwconncomp(logical(DAPI_mask));
            DAPI_mask = labelmatrix(cc);
                                   
            clear BW BW_DAPI BW_Fluo DAPI_fix DAPI_temp Fluo_fix new_cell new_nucl norm
            saved_flag=false;
        elseif button == 'c'  %% Cut - Cut nuclei using watershed algorithm
            
            Cell_num = DAPI_MaxP(uint16(round(y)),uint16(round(x)));
            
            if Cell_num ~= 0

                [i, j, k] = ind2sub([exp.image_size,exp.image_size,size(DAPI_mask,3)], cc.PixelIdxList{Cell_num});
                Cell = DAPI_mask(min(i):max(i), min(j):max(j), min(k):max(k)) == Cell_num;
                MaxP_cell_unpadded = max(Cell,[],3);

                %Pad MaxP with zeros
                MaxP_cell = zeros(size(MaxP_cell_unpadded,1)+2, size(MaxP_cell_unpadded,2)+2);
                MaxP_cell(2:(end-1),2:(end-1)) = MaxP_cell_unpadded;

                %%  Select centers of cells
                Centers_2D = zeros(size(MaxP_cell));
                %Centers_3D = zeros(size(Cell));
                figure(4)
                imagesc(MaxP_cell)
                axis equal
                [x,y,button] = ginput;
                x = uint8(round(x));
                y = uint8(round(y));


                %% Solve Laplace's equation within
                result    = MaxP_cell;
                idx       = find(MaxP_cell);
                grid      = zeros(size(MaxP_cell));
                grid(idx) = 1:length(idx);
                [M,N]     = size(grid);

                % Use the boundary pixels to form the right side of the linear system.
                boundaryCond = zeros(size(MaxP_cell));
                boundaryCond(~MaxP_cell) = -1;

                % Add in centers
                for pts = 1:length(x)
                    boundaryCond(y(pts),x(pts)) = -10;
                end

                rightside = zeros(M,N);
                rightside(2:(M-1),2:(N-1)) = boundaryCond(1:(M-2),2:(N-1)) + ...
                    boundaryCond(3:M,2:(N-1)) + boundaryCond(2:(M-1),1:(N-2)) + ...
                    boundaryCond(2:(M-1),3:N);
                rightside = rightside(idx);

                % Form the sparse D matrix from the numbered nodes of the grid matrix.
                % This part is borrowed from toolbox/matlab/demos/delsq.m.
                % Connect interior points to themselves with 4's.
                i_d = grid(idx);
                j_d = grid(idx);
                s = 4*ones(size(idx));

                % for k = north, east, south, west
                for k_d = [-1 M 1 -M]
                    % Possible neighbors in the k-th direction
                    Q = grid(idx+k_d);
                    % Index of points with interior neighbors
                    q = find(Q);
                    % Connect interior points to neighbors with -1's.
                    i_d = [i_d; grid(idx(q))]; %#ok<AGROW>
                    j_d = [j_d; Q(q)];
                    s = [s; -ones(length(q),1)];
                end
                D = sparse(i_d,j_d,s);

                % Solve the linear system.
                x = D \ rightside;
                result(idx) = x;

                result(~MaxP_cell) = -Inf;
                L = watershed(result);
                BG_label = L(1,1);
                L(L==BG_label) = 0;

                L_unpadded = L(2:(end-1),2:(end-1));

                for k_z = 1:size(Cell,3)
                    Cell_z = Cell(:,:,k_z);
                    Cell_z(L_unpadded==0) = 0;
                    Cell(:,:,k_z) = Cell_z;
                end

                DAPI_mask(cc.PixelIdxList{Cell_num}) = 0;

                DAPI_mask(min(i):max(i), min(j):max(j), min(k):max(k)) = DAPI_mask(min(i):max(i), min(j):max(j), min(k):max(k)) + uint8(Cell);

                cc = bwconncomp(logical(DAPI_mask));
                DAPI_mask = labelmatrix(cc);
                DAPI_MaxP = max(DAPI_mask, [], 3);
                [ Ann_mask ] = CreateAnnulus(DAPI_MaxP, LcFull);

                close(4)

                clear Cell Cell_num Cell_z Centers_2D D L L_unpadded MaxP_cell MaxP_cell_unpadded boundaryCond grid idx x s result rightside
                clear Q q x y button M N i_d j_d i j k k_d k_z pts
            
            end
            saved_flag=false;
        elseif button == 'd' %% Delete cell
            DAPI_Cell_num = uint16(DAPI_MaxP(uint16(round(y)),uint16(round(x))));
            Cell_num = uint16(LcFull(uint16(round(y)),uint16(round(x))));
            
            if DAPI_Cell_num ~=0
                DAPI_mask(DAPI_mask==DAPI_Cell_num) = 0;
            end     
            if Cell_num ~=0
                LcFull(LcFull==Cell_num) = 0;
            end 
                
            cc = bwconncomp(logical(DAPI_mask));
            DAPI_mask = labelmatrix(cc);
            DAPI_MaxP = max(DAPI_mask, [], 3);
            
            [ DAPI_mask, DAPI_MaxP] = LabelMatchNuc2Cell( DAPI_mask,DAPI_MaxP,LcFull );
            [ Ann_mask ] = CreateAnnulus(DAPI_MaxP, LcFull);
            
            saved_flag=false;
        elseif button == 'x'  %% Cut - Cut out cell and nuclei
            
            figure(4)
            BW = roipoly(RGB_DIC_outline);
            close(4)
            
            for i2=1:size(DAPI_mask,3)
                plane = DAPI_mask(:,:,i2);
                plane(BW) = 0;
                DAPI_mask(:,:,i2) = plane;
                Cell_mask(BW) = 0;
            end
            
            cc = bwconncomp(logical(DAPI_mask));
            DAPI_mask = labelmatrix(cc);
            DAPI_MaxP = max(DAPI_mask, [], 3);
            
            [ DAPI_mask, DAPI_MaxP] = LabelMatchNuc2Cell( DAPI_mask,DAPI_MaxP,LcFull );
            
            clear plane BW
            saved_flag=false;
        elseif button == 'f'  %% Add to cell boundary
            
            figure(4)
            BW = roipoly(RGB_DIC_outline);
            close(4)
            
%             LcFull = bwlabel(Cell_mask);
            
            cells_selected = unique(nonzeros(LcFull(BW)));
            
            if size(cells_selected,1)==1
                LcFull(BW) = cells_selected(1);
                Cell_mask(BW) = 1;                
            else
                Cell_mask(BW) = 1;
                LcFull(BW) = max(LcFull(:))+1;
            end
            
            cc = bwconncomp(logical(DAPI_mask));
            DAPI_mask = labelmatrix(cc);
            DAPI_MaxP = max(DAPI_mask, [], 3);
            [ DAPI_mask, DAPI_MaxP] = LabelMatchNuc2Cell( DAPI_mask,DAPI_MaxP,LcFull );
            
            [ Ann_mask ] = CreateAnnulus(DAPI_MaxP, LcFull);
            
            clear BW cells_selected 
            saved_flag=false;
        elseif button == 'w'  %% Use DIC to define cell boundary - looks for white cell boundary.
                
            
            DIC_temp = imread([exp.image_dir exp.name_root num2str(image_id,'%03d') 'z' num2str(6,['%0' num2str(exp.dig_z) 'd']) 'c1.tif']);

            DIC_temp = imadjust(DIC_temp);
            
            %DIC_temp = imclose(DIC_temp, strel('disk', 10));
            
            DIC_temp = imfilter(DIC_temp, fspecial('log', 10));
            
            DIC_temp(1:10,:) = 0;
            DIC_temp(:,1:10) = 0;
            DIC_temp(1014:1024,:) = 0;
            DIC_temp(:,1014:1024) = 0;
            
            DIC_temp = im2bw(DIC_temp, graythresh(DIC_temp));
            
            DIC_temp = imclose(DIC_temp, strel('disk', 10));
            
            DIC_temp = imfill(DIC_temp, 'holes');
            
            Cell_mask = DIC_temp | Cell_mask;
            
            clear DIC_temp
                         
            
            saved_flag=false;
        elseif button == 'g'  %% Get rid over overlapping cells, discarding smaller cells
            
            %Find the order of cells, smallest to largest
            max_cell = single(max(DAPI_mask(:)));
            pix = single(DAPI_mask(:));
            n_pix_cell = hist(pix,0:max_cell);
            [sort_cell_size, i_cell] = sort(n_pix_cell, 'ascend');
            i_cell = i_cell - 1; % To convert from position of hist to pix id.
            
            %Find the x-y indicies for all cells.
            for cell_num = i_cell(i_cell~=0) 
                [x,y,z] = ind2sub(size(DAPI_mask), find(DAPI_mask==cell_num));
                Ind{cell_num} = sub2ind([size(DAPI_mask,1),size(DAPI_mask,2)], x, y)';
            end
               
            %For each, check if there's overlap with larger cells
            for i__ = 1:(numel(i_cell)-2)
                pix_of_larger_cells = [Ind{i_cell((i__+1):(end-1))}];
                overlap = intersect(Ind{i_cell(i__)}',pix_of_larger_cells');
                if ~isempty(overlap)
                    DAPI_mask(DAPI_mask==i_cell(i__)) = 0;
                end
            end
            
            cc = bwconncomp(logical(DAPI_mask));
            DAPI_mask = labelmatrix(cc);
            DAPI_MaxP = max(DAPI_mask, [], 3);
            
            clear pix_of_larger_cells overlap i__ i_cell sort_cell_size n_pix_cell pix max_cell Ind x y z cell_num
            saved_flag=false;
        elseif button == 'r'  %% Recalculate threshold for cell boundary
            
            slices = size(DAPI_mask,3);
            
            Fluo_temp = zeros(exp.image_size,exp.image_size,slices);
            x_Fluo = 1E2:1E2:1E5;
            
            for i2 = 1:slices;
                Fluo_temp(:,:,i2) = double(imread([exp.image_dir exp.name_root num2str(image_id,'%03d') 'z' num2str(i2,['%0' num2str(exp.dig_z) 'd']) 'c' num2str(exp.Fluo_channel) '.tif']));                
            end
            
            Fluo_image_maxp = max(Fluo_temp,[],3);
            
            n_Fluo = hist(Fluo_image_maxp(:),x_Fluo);
            f_Fluo = n_Fluo./numel(Fluo_image_maxp(:));
            
            figure(4)
            plot(x_Fluo, f_Fluo);
            set(gca, 'YScale', 'log')
            [x,y]=ginput(1);
            close(4);
            disp([sprintf('\n') 'New threshold = ' num2str(x)])
            
            Cell_mask = zeros(exp.image_size,exp.image_size);
            Cell_mask(Fluo_image_maxp>=x) = 1;
            Cell_mask = imfill(Cell_mask,'holes');
            Cell_mask = imopen(Cell_mask,strel('disk',5));
            
            Cell_mask(DAPI_MaxP~=0) = 1;
            
            cc = bwconncomp(logical(DAPI_mask));
            DAPI_mask = labelmatrix(cc);
            DAPI_MaxP = max(DAPI_mask, [], 3);
            
            clear x y Fluo_image_maxp Fluo_temp x_Fluo n_Fluo f_Fluo
                        
            % AF added 20151024
            saved_flag=false;
        elseif button == 'j'  %% Join two cells together                                    
        %%  Select cells to join                            
            figure(4); imagesc(LcFull)
%                 axis equal
            inputloop=1;
%                 disk_size=2;
            while inputloop
                [x,y,button] = ginput;                                
                if length(x) == 2                                                                                          
                    inputloop=0;                    
                else
                    disp('Please connect only two cells')
                end
            end
            clear inputloop

            x = uint16(round(x));
            y = uint16(round(y));

            Cell1=zeros(size(LcFull));
            Cell2=zeros(size(LcFull));

            Cell1(LcFull==LcFull(y(1),x(1)))=1;
            Cell2(LcFull==LcFull(y(2),x(2)))=1;

            Cell1_prm=bwperim(Cell1); 
            Cell2_prm=bwperim(Cell2);

            a=find(Cell1_prm);
            b=find(Cell2_prm);
            % Determine pairwise distance for every permimeter point
            for i=1:numel(a)
                [x1,y1]=ind2sub(size(Cell1),a(i));                            
                for j=1:numel(b)
                    [x2,y2]=ind2sub(size(Cell1),b(j));
                    d2(j)=sqrt((x2-x1)^2+(y2-y1)^2);
                end
                d1(i)=min(d2);
                d_list(i,1)=a(i);
                d_list(i,2)=b(find(d2==d1(i),1));
                clear d2
            end
            clear a b i j x1 x2 y1 y2

            % draw lines between points of the cells that are close together
            a=find(d1<min(d1)*1.5+1);
            Cell_join=Cell1_prm | Cell2_prm;
            for i=1:numel(a)
                [y1,x1]=ind2sub(size(Cell1),d_list(a(i),1));
                [y2,x2]=ind2sub(size(Cell1),d_list(a(i),2));

                % distances according to both axes
                xn = abs(x2-x1);
                yn = abs(y2-y1);

                % interpolate against axis with greater distance between points;
                % this guarantees statement in the under the first point!
                if (xn > yn)
                    xc = x1 : sign(x2-x1) : x2;
                    yc = round( interp1([x1 x2], [y1 y2], xc, 'linear') );
                else
                    yc = y1 : sign(y2-y1) : y2;
                    xc = round( interp1([y1 y2], [x1 x2], yc, 'linear') );
                end

                % 2-D indexes of line are saved in (xc, yc), and
                % 1-D indexes are calculated here:
                ind = sub2ind( size(Cell1), yc, xc );

                % draw line on the image
                Cell_join(ind) = 1;
            end
            Cell_join=imfill(Cell_join,'holes');

            Cell_mask = LcFull | logical(Cell_join); 
            LcFull = bwlabel(Cell_mask);            
            %  Match label number of nuclei to label number of cells
            [ DAPI_mask, DAPI_MaxP] = LabelMatchNuc2Cell( DAPI_mask,DAPI_MaxP,LcFull );
            % update Annuli
            [ Ann_mask ] = CreateAnnulus(DAPI_MaxP, LcFull);
            
            close(4)
            saved_flag=false;
        elseif button == 'k'  %% Cut Cell mask
            figure(4)
            
            close(4)
            
            saved_flag=false;
        elseif button == 'l'  %% Reload from file
            load([exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_2.mat'], 'DAPI_mask', 'LcFull', 'Ann_mask') 
            saved_flag=true;
        elseif button == 's' %% Save masks            
            disp(['Saving masks to ' exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_2.mat\n'])
            save([exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_2.mat'], 'DAPI_mask', 'LcFull', 'Ann_mask')   
            saved_flag=true;
        end            
    end
    
    
    %%  Match label number of nuclei to label number of cells
    [ DAPI_mask, DAPI_MaxP] = LabelMatchNuc2Cell( DAPI_mask,DAPI_MaxP,LcFull );
    [ Ann_mask ] = LabelMatchAnn2Cell( Ann_mask,LcFull );
    
    
    
    
%     figure(1)
%     patch(isosurface(DAPI_mask,.5),'FaceColor','blue','EdgeColor','none');
%     view(3); axis equal
%     camlight; lighting phong
% 
%     figure(2)
%     patch(isosurface(InNu_mask,.5),'FaceColor','blue','EdgeColor','none');
%     view(3); axis equal
%     camlight; lighting phong
    
    if ~saved_flag
        if strcmp(save_button,'Yes')
            disp(['Saving masks to ' exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_2.mat\n'])
            save([exp.root_dir 'Data Analysis\Segmentation\' exp.name_root 'seg' num2str(image_id, '%03d') '_2.mat'], 'DAPI_mask',  'LcFull', 'Ann_mask') 
        end      
    end
    
    clear Cell Centers D DAPI L mask
    close all
end




%% Unused Old Code
% %% Create inner nucleus mask by shrinking manually fixed nucleus mask
%     InNu_mask = zeros(size(DAPI_mask));
%     DAPI_cc = bwconncomp(DAPI_mask);
%     DAPI_rp = regionprops(DAPI_cc, 'SubarrayIdx', 'PixelIdxList', 'Centroid');
%        
%     for i2 = 1:DAPI_cc.NumObjects
%         
%         xysize = 25;
%         zsize = 4;
%         
%         center = double(round(DAPI_rp(i2).Centroid));
%         
%         center(3) = min([size(DAPI_mask,3)-5 center(3)]);
%         center(3) = max([5 center(3)]);
%         
%         if ((exp.image_size/2)-abs(((exp.image_size/2)+0.5)-center(1))+.5)-1<xysize; xysize = double((exp.image_size/2)-abs(((exp.image_size/2)+0.5)-center(1))+.5)-1; end
%         
%         if ((exp.image_size/2)-abs(((exp.image_size/2)+0.5)-center(2))+.5)-1<xysize; xysize = double((exp.image_size/2)-abs(((exp.image_size/2)+0.5)-center(2))+.5)-1; end
%         
%         if ( 10-abs( 10.5-center(3))+.5)-1<zsize;   zsize = double( 10-abs( 10.5-center(3))+.5)-1; end
% 
%         [ellipse1] = makeellipse(xysize,zsize);
%     
%         InNu_mask((center(2)-xysize):(center(2)+xysize),(center(1)-xysize):(center(1)+xysize),(center(3)-zsize):(center(3)+zsize)) = i2*ellipse1;
%         
%     end
%     InNu_mask = logical(InNu_mask).*logical(DAPI_mask);
%     clear  DAPI_rp ellipse1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  %% Cell outlines by thresholding
%  close all
%  figure(1); imshow(mat2gray(Disp))
%  %%
%     clear temp temp2 temp3 temp4 DIC_temp thesh
%     DIC_temp = imread([exp.image_dir exp.name_root num2str(image,'%03d') 'z' num2str(6,['%0' num2str(exp.dig_z) 'd']) 'c1.tif']);
% 
%     DIC_temp = imadjust(DIC_temp);
%     
% %%
%     %DIC_temp = imclose(DIC_temp, strel('disk', 10));
% 
% %     DIC_temp = imfilter(DIC_temp, fspecial('log', 10));
% %     
% %     % Crop borders
% %     DIC_temp(1:10,:) = 0;
% %     DIC_temp(:,1:10) = 0;
% %     DIC_temp(1014:1024,:) = 0;
% %     DIC_temp(:,1014:1024) = 0;
% % %%
%     DIC_thresh = im2bw(DIC_temp, graythresh(DIC_temp));
%     DIC_thresh = imerode(DIC_thresh,strel('disk',5));
%     DIC_thresh = imdilate(DIC_thresh,strel('disk',2));
% %     DIC_thresh = imclose(DIC_thresh,strel('disk',10));
%     figure(2); imshow(DIC_thresh)
%     %%
% %     [temp,thresh] = edge(DIC_temp,'Canny');
% %     temp = edge(DIC_temp,'Canny',[thresh(1) 4*thresh(2)]);
% %     temp = imdilate(temp,strel('disk',5));
% %     figure(3); imshow(temp)
%     %
%     temp2 = imclose(DIC_thresh,strel('disk',30));
% %     temp2 = ~temp2;
%     %
% %     temp2=imclose(temp2,strel('disk', 25));
%     temp2=imfill(temp2,'holes');
%     figure(4); imshow(temp2);
%     %%
% %     DIC_thresh = imerode(DIC_thresh,strel('disk',5));
% %     DIC_thresh = imdilate(DIC_thresh,strel('disk',2));
% %     DIC_thresh = imclose(DIC_thresh,strel('disk',15));
%     
%     temp3= temp2 & ~DIC_thresh;
%     temp3=imopen(temp3,strel('disk',8));
%     temp3=imfill(temp3,'holes');
% %     temp3= temp3 & ~temp;
%     figure(5); imshow(temp3);
%     %%
% %     temp4= temp2 & ~temp & ~DIC_thresh;
% %     
% %     temp4=imopen(temp4,strel('disk',5));
% %     temp4=imclose(temp4,strel('disk',10));
% %     
% %     temp4=imfill(temp4,'holes');
% % %     temp3= temp3 & ~temp;
% %     figure(6); imshow(temp4);
%     
%     %%
%     figure(7); imshowpair(mat2gray(Disp),temp3);
% %     figure(8); imshowpair(mat2gray(Disp),temp4);
%     
%     %%
% %     DIC_temp = edge(DIC_temp);
% 
%     DIC_temp = imclose(temp, strel('disk', 15));
%     
% 
%     DIC_temp = imfill(DIC_temp, 'holes');
%     
%     DAPI_MaxP = max(DAPI_mask, [], 3);
%     Cell_mask = DIC_temp | RFP_MaxP;
%     Cell_mask = DAPI_MaxP | Cell_mask;
% 
%     clear DIC_temp         
