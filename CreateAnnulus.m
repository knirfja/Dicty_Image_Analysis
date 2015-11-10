function [ Ann_mask ] = CreateAnnulus(DAPI_MaxP, LcFull)
%CreateAnnulus Summary of this function goes here
%   Detailed explanation goes here


% Create annulus
Ann_mask = repmat(uint8(0),size(DAPI_MaxP,1),size(DAPI_MaxP,2));
for iter=nonzeros(unique(DAPI_MaxP(:)))'
%         disp(['\nCreating annulus for ' num2str(iter)])
    Ann_pix = Ann_mask==iter;
    Nucl_pix = DAPI_MaxP==iter;

    Ann_temp = repmat(uint8(0),size(Ann_mask,1),size(Ann_mask,2));
    Ann_temp(Nucl_pix) = iter;
    Ann_temp = imdilate(Ann_temp, strel('disk',10));
    Ann_temp = imfill(Ann_temp,'holes');
    Ann_temp(Nucl_pix)=0; % subtract nucl
    Ann_temp(LcFull~=iter)=0; % subtract pix outside of cell        

    if sum(LcFull(:)==iter) <= 2*sum(Nucl_pix(:))
        % If cell is too small, use entire cell
        Ann_temp(LcFull==iter)=iter;
        Ann_temp(logical(DAPI_MaxP(:))) = 0; % subtract nucl
        Ann_temp(LcFull~=iter)=0; % subtract pix outside of cell
    else
        % dilate until at least the area of nuc        
        while sum(Ann_pix(:)) <= sum(Nucl_pix(:))          
            Ann_temp = imdilate(Ann_temp,strel('disk',10));
            Ann_temp(logical(DAPI_MaxP(:))) = 0; % subtract nucl
            Ann_temp(LcFull~=iter)=0; % subtract pix outside of cell
            Ann_pix = Ann_temp==iter;            
        end
    end
    Ann_mask(logical(Ann_temp))=iter;        
end
clear iter Ann_temp Ann_pix Nucl_pix  

end

