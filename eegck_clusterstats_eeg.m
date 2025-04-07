function [PosClus,NegClus] = eegck_clusterstats_eeg(cfg,X,XR)

% Compute significant clusters in array X, based on the shuffled data XR
% this version assumes that the data have dimension EEG_CHan X time (or freq)
% 
% function [clusters] = eegck_clusterstats_eeg(cfg,X,XR)
%
%
% cfg.critvaltype ='par'; 'prctile' % type of threshold to apply. Usual 'par'
% cfg.critval = 2 ; critical cutoff value for cluster members if parametric
% cfg.critval = [2.5, 97.5]  for prctiles  
% cfg.conn = 8; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
% cfg.clusterstatistic = 'max' 'maxsize' 'maxsum'
% cfg.minsize = 2; % minimal cluster size
% cfg.pval = 0.05; % threshold to select signifciant clusters
% cfg.df = degrees of freedom. Only needed for effect size.
% cfg.neighbours  -> FT neighbourhood structure
%
%
% output: 
% .p  p-values for each cluster
% .stat  cluster stats for each cluster, e.g. summed cluster statistice for maxsum
% .m_Effect  maximal (poscluster) minimal (negcluster) effect within a cluster
% maskSig mask with only signifciant clusters


nel = size(X);

% neighbourbood matrix to index
labels  = {cfg.neighbours.label};
NBMat = zeros(nel(1),nel(1));
for e=1:nel(1)
  nb = cfg.neighbours(e).neighblabel;
  for b=1:length(nb)
    j = find(strcmp(cfg.neighbours(e).neighblabel{b},labels));
    if ~isempty(j)
      NBMat(e,j) =1;
    end
  end
end

    
 ft_hastoolbox('spm8',1);
%-----------------------------------------------------------
% switch stats
switch cfg.critvaltype
  case {'par','Par'}
    % keep as is, double for pos / neg
    cfg.critval = [-cfg.critval cfg.critval];
    
  case {'Prctile','prctile'}
    cfg.critval(1) = prctile(XR(:),cfg.critval(1));
    cfg.critval(2) = prctile(XR(:),cfg.critval(2));
end

%----------------------------------------------------------------
% positive true clusters
tmp = double( (X>=cfg.critval(2)));
% connect neighbouring channels
% posclusobs = findclusterck(reshape(postailrnd(:,i), [cfg.dim,1]),NBMat,cfg.minsize);
posclusobs = findclusterck(tmp,NBMat,cfg.minsize);
% collect stats
Nobspos = max(posclusobs(:));
statPO = zeros(Nobspos,1);

for j = 1:Nobspos
  if strcmp(cfg.clusterstatistic, 'max'),
    statPO(j) = max(X(posclusobs==j));
  elseif strcmp(cfg.clusterstatistic, 'maxsize'),
    statPO(j) = length(find(posclusobs==j));
  elseif strcmp(cfg.clusterstatistic, 'maxsum'),
    statPO(j) = nansum(X(posclusobs==j));
  end
  EffectPO(j) = max(X(posclusobs==j));
end

%----------------------------------------------------------------
% negative true clusters
tmp = double( (X<=cfg.critval(1)));
negclusobs = findclusterck(tmp,NBMat,cfg.minsize);
Nobsneg = max(negclusobs(:));
statNO = zeros(Nobsneg,1);

for j = 1:Nobsneg
  if strcmp(cfg.clusterstatistic, 'max'),
    statNO(j) = min(X(negclusobs==j));
  elseif strcmp(cfg.clusterstatistic, 'maxsize'),
    statNO(j) = -length(find(negclusobs==j));
  elseif strcmp(cfg.clusterstatistic, 'maxsum'),
    statNO(j) = nansum(X(negclusobs==j));
  end
  EffectNO(j) = min(X(negclusobs==j));
end



%----------------------------------------------------------------------------------------------------------
% process randomized data

Nrand = size(XR,3);
StatRand = zeros(Nrand,2);
for i=1:Nrand
  % positive clusters
  tmp = double( (XR(:,:,i)>=cfg.critval(2)));
  posclusR = findclusterck(tmp,NBMat,cfg.minsize);
  
  % negative clusters
  tmp = double( (XR(:,:,i)<=cfg.critval(1)));
  [negclusR, negnum] = bwlabeln(tmp);
  
  % cluster-statistics
  tmp = squeeze(XR(:,:,i));
  Nobspos = max(posclusR(:));
  statPR = zeros(3,Nobspos);
  % collect the different statistics
  for j = 1:Nobspos
    statPR(1,j) = max(tmp(posclusR==j));
    statPR(2,j) = length(find(posclusR==j));
    statPR(3,j) = nansum(tmp(posclusR==j));
  end
  
  
  Nobsneg = max(negclusR(:));
  statNR = zeros(3,Nobsneg);
  for j = 1:Nobsneg
    statNR(1,j) = min(tmp(negclusR==j));
    statNR(2,j) = -length(find(negclusR==j));
    statNR(3,j) = nansum(tmp(negclusR==j));
  end
  
  % select stats
  if strcmp(cfg.clusterstatistic, 'max'),
    statNR = statNR(1,:);
    statPR = statPR(1,:);
  elseif strcmp(cfg.clusterstatistic, 'maxsize'),
    statNR = statNR(2,:);
    statPR = statPR(2,:);
  elseif strcmp(cfg.clusterstatistic, 'maxsum'),
    statNR = statNR(3,:);
    statPR = statPR(3,:);
  end
  % max-stats
  pout = 0;
  nout =0;
  if ~isempty(statPR)
    pout = max(statPR);
  end
  if ~isempty(statNR)
    nout = min(statNR);
  end
  StatRand(i,:) = [pout nout]';
end


%----------------------------------------------------------------------------------------------------------
% compare actual and random distributions
PosClus=[];
Npos = length(statPO);
for i=1:Npos
  PosClus.p(i) = nansum(StatRand(:,1)>statPO(i))/Nrand;
  PosClus.stat(i) = statPO(i);
  PosClus.mask = double(posclusobs);
  if PosClus.p(i)==0
    PosClus.p(i)=1/Nrand;
  end
  PosClus.Effect(i) = EffectPO(i);
end

NegClus=[];
Nneg = length(statNO);
for i=1:Nneg
  NegClus.p(i) = nansum(StatRand(:,2)<statNO(i))/Nrand;
  NegClus.stat(i) = statNO(i);
  NegClus.mask = double(negclusobs);
  if NegClus.p(i)==0
    NegClus.p(i)=1/Nrand;
  end
  NegClus.Effect(i) = EffectNO(i);
end

%----------------------------------------------------------------------------------------------------------
% remove nonsig clusters from mask

if ~isempty(PosClus)
  PosClus.maskSig = PosClus.mask;
  J = find(PosClus.p>cfg.pval);
  for k=1:length(J)
    PosClus.maskSig(find(PosClus.maskSig(:)==J(k))) = 0;
  end
end

if ~isempty(NegClus)
  
  J = find(NegClus.p>cfg.pval);
  
  NegClus.maskSig = NegClus.mask;
  for k=1:length(J)
    NegClus.maskSig(find(NegClus.maskSig(:)==J(k))) = 0;
  end
end





return;



%%

function [cluster, total] = findclusterck(onoff, spatdimneighbstructmat, varargin)

%  findclusterck(tmp,NBMat,cfg.minsize);

% FINDCLUSTER returns all connected clusters in a 3 dimensional matrix
% with a connectivity of 6.
%
% Use as
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, minnbchan)
% or ar
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, spatdimneighbselmat, minnbchan)
% where 
%   onoff                   is a 3D boolean matrix with size N1xN2xN3
%   spatdimneighbstructmat  defines the neighbouring channels/combinations, see below
%   minnbchan               the minimum number of neighbouring channels/combinations 
%   spatdimneighbselmat     is a special neighbourhood matrix that is used for selecting
%                           channels/combinations on the basis of the minnbchan criterium
%
% The neighbourhood structure for the first dimension is specified using 
% spatdimneighbstructmat, which is a 2D (N1xN1) matrix. Each row and each column corresponds
% to a channel (combination) along the first dimension and along that row/column, elements
% with "1" define the neighbouring channel(s) (combinations). The first dimension of
% onoff should correspond to the channel(s) (combinations).
%
% See also SPM_BWLABEL (spm toolbox) 

% Copyright (C) 2004, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

spatdimlength = size(onoff, 1);
nfreq = size(onoff, 2);
ntime = size(onoff, 3);

if length(size(spatdimneighbstructmat))~=2 || ~all(size(spatdimneighbstructmat)==spatdimlength)
  ft_error('invalid dimension of spatdimneighbstructmat');
end

minnbchan=0;
if length(varargin)==1
    minnbchan=varargin{1};
end
if length(varargin)==2
    spatdimneighbselmat=varargin{1};
    minnbchan=varargin{2};
end

if minnbchan>0
    % For every (time,frequency)-element, it is calculated how many significant
    % neighbours this channel has. If a significant channel has less than minnbchan
    % significant neighbours, then this channel is removed from onoff.
    
    if length(varargin)==1
        selectmat = single(spatdimneighbstructmat | spatdimneighbstructmat');
    end
    if length(varargin)==2
        selectmat = single(spatdimneighbselmat | spatdimneighbselmat');
    end
    nremoved=1;
    while nremoved>0
      nsigneighb=reshape(selectmat*reshape(single(onoff),[spatdimlength (nfreq*ntime)]),[spatdimlength nfreq ntime]);
      remove=(onoff.*nsigneighb)<minnbchan;
      nremoved=length(find(remove.*onoff));
      onoff(remove)=0;
    end
end

% for each channel (combination), find the connected time-frequency clusters
labelmat = zeros(size(onoff));
total = 0;
if nfreq*ntime>1
  for spatdimlev=1:spatdimlength
    [labelmat(spatdimlev, :, :), num] = spm_bwlabel(double(reshape(onoff(spatdimlev, :, :), nfreq, ntime)), 6); % the previous code contained a '4' for input
    labelmat(spatdimlev, :, :) = labelmat(spatdimlev, :, :) + (labelmat(spatdimlev, :, :)~=0)*total;
    total = total + num;
  end
else
  labelmat(onoff>0) = 1:sum(onoff(:));
  total = sum(onoff(:));
end
% combine the time and frequency dimension for simplicity
labelmat = reshape(labelmat, spatdimlength, nfreq*ntime);

% combine clusters that are connected in neighbouring channel(s)
% (combinations). Convert inputs to uint32 as that is required by the mex
% file (and the values will be positive integers anyway).
cluster = combineClusters(uint32(labelmat), logical(spatdimneighbstructmat), uint32(total));

% reshape the output to the original format of the data
cluster = reshape(cluster, spatdimlength, nfreq, ntime);

% update the total number
total = numel(unique(cluster(:)))-1; % the value of 0 does not count
return;