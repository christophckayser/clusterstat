function [PosClus,NegClus] = eegck_clusterstats(cfg,X,XR)

% Compute significant clusters in array X, based on the shuffled data XR
% Assumes that the array reflects continous and connected data, e.g. time, frequency or source space
% 
% function [PosClus,NegClus] = eegck_clusterstats(cfg,X,XR)
%
% X can be 1D, 2D or 3D; XR the same + randomization dimension (last)
%
%
% cfg.critvaltype ='par'; 'prctile' % type of threshold to apply. Usual 'par'
% cfg.critval = 2 ; critical cutoff value for cluster members if parametric
% cfg.critval = [2.5, 97.5]  for prctiles  
% 
% cfg.conn = 6; % connectivity criterion (for 2D 4 or 8 for 3D 6,18,26) 
% cfg.clusterstatistic = 'max' 'maxsize' 'maxsum'
% cfg.minsize = 2; % minimal cluster size
% cfg.pval = 0.05; % threshold to select signifciant clusters
% cfg.df = degrees of freedom. Only needed for effect size.
%
% output: 
% .p  p-values for each cluster
% .stat  cluster stats for each cluster, e.g. summed cluster statistice for maxsum
% .m_Effect  maximal (poscluster) minimal (negcluster) effect within a cluster
% mask  mask with cluster IDs, both signifcant and not signif
% maskSig mask with only signifciant clusters


% determine dimensions
if min(size(X))==1
  dim = 1;
  X = X(:);
else
  dim = ndims(X);
end

 ft_hastoolbox('spm8',1);
%-----------------------------------------------------------
% process actual data

switch cfg.critvaltype
  case {'par','Par'}
    % keep as is, double for pos / neg
    cfg.critval = [-cfg.critval cfg.critval];
    
  case {'Prctile','prctile'}
    cfg.critval(1) = prctile(XR(:),cfg.critval(1));
    cfg.critval(2) = prctile(XR(:),cfg.critval(2));
end

% pos
tmp = double( (X>cfg.critval(2)));
if dim>1
  [posclusobs,~] = spm_bwlabel(tmp,cfg.conn);
else
   [posclusobs] = bwlabeln(tmp,cfg.conn);
end
% collect stats
Nobspos = max(posclusobs(:));
statPO = zeros(5,Nobspos);

for j = 1:Nobspos
  % max
  statPO(1,j) = max(X(posclusobs==j));
  % size
  statPO(2,j) = length(find(posclusobs==j));
  % sum
  statPO(3,j) = nansum(X(posclusobs==j));
  statPO(4,j) = j;
  % max-value
  statPO(5,j) = max(X(posclusobs==j));

end

% neg
tmp = double( (X<cfg.critval(1)));
if dim>1
  [negclusobs,~] = spm_bwlabel(tmp,cfg.conn);
else
   [negclusobs] = bwlabeln(tmp,cfg.conn);
end
Nobsneg = max(negclusobs(:));
statNO = zeros(5,Nobsneg);
for j = 1:Nobsneg
  statNO(1,j) = min(X(negclusobs==j));
  statNO(2,j) = -length(find(negclusobs==j));
  statNO(3,j) = nansum(X(negclusobs==j));
  statNO(4,j) = j;
  statNO(5,j) = min(X(negclusobs==j));
  
end
% threshold sizes
if cfg.minsize,
  J = find( statPO(2,:)>= cfg.minsize);
  statPO = statPO(:,J);
  posclusobs2 = posclusobs*0;
  for l=1:length(J)
    posclusobs2(find(posclusobs(:)==statPO(4,l))) = l;
  end
  
  J = find( statNO(2,:)<= -cfg.minsize);
  statNO = statNO(:,J);
  negclusobs2 = negclusobs*0;
  for l=1:length(J)
    negclusobs2(find(negclusobs(:)==statNO(4,l))) = l;
  end
end

% max / min values
if length(statPO)
  EffectPos = statPO(5,:);
end
if length(statNO)
  EffectNeg = statNO(5,:);
end

% select stats
if strcmp(cfg.clusterstatistic, 'max'),
  statNO = statNO(1,:);
  statPO = statPO(1,:);
elseif strcmp(cfg.clusterstatistic, 'maxsize'),
  statNO = statNO(2,:);
  statPO = statPO(2,:);
elseif strcmp(cfg.clusterstatistic, 'maxsum'),
  statNO = statNO(3,:);
  statPO = statPO(3,:);
end

    
%----------------------------------------------------------------------------------------------------------
% process randomized data

ndimxr = size(XR);
Nrand = ndimxr(end);
StatRand = zeros(Nrand,2);
for i=1:Nrand
  if dim==1
    tmp = double( (XR(:,i)>cfg.critval(2)));
    [posclusR, posnum] = bwlabeln(tmp,cfg.conn);
    % neg
    tmp = double( (XR(:,i)<cfg.critval(1)));
    [negclusR, negnum] = bwlabeln(tmp,cfg.conn);
    tmp = sq(XR(:,i));
  elseif dim==2
    tmp = double( (XR(:,:,i)>cfg.critval(2)));
    [posclusR, posnum] = spm_bwlabel(tmp,cfg.conn);
    % neg
    tmp = double( (XR(:,:,i)<cfg.critval(1)));
    [negclusR, negnum] = spm_bwlabel(tmp,cfg.conn);
    tmp = sq(XR(:,:,i));
  elseif dim==3
    tmp = double( (XR(:,:,:,i)>cfg.critval(2)));
    [posclusR, posnum] = spm_bwlabel(tmp,cfg.conn);
    % neg
    tmp = double( (XR(:,:,:,i)<cfg.critval(1)));
    [negclusR, negnum] = spm_bwlabel(tmp,cfg.conn);
    tmp = sq(XR(:,:,:,i));
  end
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
  
  % threshold sizes
  if cfg.minsize,
    J = find( statPR(2,:)>= cfg.minsize);
    statPR = statPR(:,J);
    J = find( statNR(2,:)<= -cfg.minsize);
    statNR = statNR(:,J);
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
  PosClus.p(i) = nansum(StatRand(:,1)>=statPO(i))/Nrand;
  PosClus.stat(i) = statPO(i);
  PosClus.mask = posclusobs2;
  PosClus.Effect(i) = EffectPos(i);
  
  if PosClus.p(i)==0
    PosClus.p(i)=1/Nrand;
  end
end

NegClus=[];
Nneg = length(statNO);
for i=1:Nneg
  NegClus.p(i) = nansum(StatRand(:,2)<=statNO(i))/Nrand;
  NegClus.stat(i) = statNO(i);
  NegClus.mask = negclusobs2;
  NegClus.Effect(i) = EffectNeg(i);
  if NegClus.p(i)==0
    NegClus.p(i)=1/Nrand;
  end
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testing for 2D





for S=1:20
  
  Act(S,:,:,1) = randn(20,40);
  Act(S,:,:,2) = randn(20,40);
end

% introduce difference
Act(:,5:6,4:10,1) = Act(:,5:6,4:10,1)+1;

Ttrue = Act(:,:,:,1)-Act(:,:,:,2);
Ttrue = sq(mean(Ttrue)./sem(Ttrue));


for boot=1:500
    Bct(S,:,:,1,boot) = randn(20,40);
    Bct(S,:,:,2,boot) = randn(20,40);
    Bct(:,5:6,4:10,1,boot) =Bct(:,5:6,4:10,1,boot)  + randn;
  
end
Tboot = sq(Bct(:,:,:,1,:)-Bct(:,:,:,2,:));
Tboot = sq(mean(Tboot)./sem(Tboot));




cfg=[];
cfg.critvaltype ='par'; %'prctile' % type of threshold to apply. Usual 'par'
cfg.critval = 2 ; %critical cutoff value for cluster members if parametric

cfg.conn = 6; % connectivity criterion (for 2D/3D)
cfg.clusterstatistic =  'maxsize' 
cfg.minsize = 4; % minimal cluster size
cfg.pval = 0.05; % threshold to select signifciant clusters
cfg.df = 20; %freedom. Only needed for effect size.

[PosClus,NegClus] = eegck_clusterstats(cfg,Ttrue,Tboot)

PosClus







