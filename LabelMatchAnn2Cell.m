function [ Ann_mask ] = LabelMatchAnn2Cell( Ann_mask,LcFull )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%  Match label number of cells to label number of nuclei
temp = zeros(size(Ann_mask));
for ann_id = int8(unique(nonzeros(Ann_mask(:))))'
    cell_id = sort(unique(nonzeros(LcFull(Ann_mask==ann_id))));
%                 LcFull(LcFull==nic_id) = -1*nuc_id;
    try
        temp(Ann_mask==ann_id) = -1*cell_id;
    catch
        ann_id
        cell_id
        error('Error, annulus mapped to multiple cells')
        
    end                                  
end
Ann_mask = -1*temp;
end

