function [mask,index] = eegck_stripclusters(PosClus,NegClus,ndim)

% function [mask,index] = eegck_stripclusters(PosClus,NegClus,ndim)
% 
% returns just the pos/ neg significant cluster-masks, e.g. to use for plotting
% ndim is the dimensions of the expected stat field, e.g. [128,1]

mask_pos = zeros(ndim);
index =[];

if ~isempty(PosClus)
  if isfield(PosClus,'maskSig')
    mask_pos = PosClus.maskSig;
  end
end

mask_neg = zeros(ndim);

if ~isempty(NegClus)
  if isfield(NegClus,'maskSig')
    mask_neg = NegClus.maskSig;
  end
end
% combine

ndims = sum(ndim>1);
if ndims==1
  mask = cat(2,mask_pos,mask_neg);
  index = cat(1,index,find(mask(:,1)));
  index = cat(1,index,find(mask(:,2)));
else
    mask = logical(mask_pos+mask_neg);
    index = sum(mask(:));
end



return;
