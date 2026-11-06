% curve_flip_test.m -- is R6's curve assignment optimal at its own plug settings?
D=load('docr_corrected_data62.mat'); L=load('deploy62.mat'); S=L.DEPLOY.R6;
PS=S.PS(:); ix0=S.ix(:); [T0,~]=eval_fixed(PS,ix0,D);
nImp=0; nInf=0; best1=T0;
for r=1:D.nR, for c=1:3, if c==ix0(r), continue; end
  ix=ix0; ix(r)=c; [Tt,Nt]=eval_fixed(PS,ix,D);
  if Nt>0, nInf=nInf+1; else, best1=min(best1,Tt); nImp=nImp+(Tt<T0-1e-9); end, end, end
PSMp=D.IfP(:)./(D.CT*PS); near=find(PSMp<6 & D.IfP(:)>0)'; K=numel(near); best=T0; bestix=ix0;
for code=0:3^K-1, ix=ix0; v=code;
  for j=1:K, ix(near(j))=mod(v,3)+1; v=floor(v/3); end
  [Tt,Nt]=eval_fixed(PS,ix,D); if Nt==0 && Tt<best, best=Tt; bestix=ix; end, end
CFT=struct('T0',T0,'nSingle',2*D.nR,'nImp',nImp,'nInf',nInf,'best1',best1,'near',near,'nComb',3^K,'bestComb',best,'nDiff',sum(bestix~=ix0));
save('curve_flip62.mat','CFT');