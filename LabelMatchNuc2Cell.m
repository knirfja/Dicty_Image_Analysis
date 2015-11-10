function [ DAPI_mask, DAPI_MaxP] = LabelMatchNuc2Cell( DAPI_mask,DAPI_MaxP,LcFull )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%  Match label number of cells to label number of nuclei
temp = zeros(size(DAPI_mask));
for nuc_id = int8(unique(nonzeros(DAPI_MaxP(:))))'
    cell_id = sort(unique(nonzeros(LcFull(DAPI_MaxP==nuc_id))));
%                 LcFull(LcFull==nic_id) = -1*nuc_id;
    try
        temp(DAPI_mask==nuc_id) = -1*cell_id;
    catch
        error('Error, nuclei mapped to multiple cells')
    end                                  
end
DAPI_mask = -1*temp;
DAPI_MaxP = max(DAPI_mask, [], 3); 
end

