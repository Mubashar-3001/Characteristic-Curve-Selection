function [Top,NV,TM] = eval_fixed(PS,ix,D)
% Exact inner solve for FIXED plug settings PS and FIXED curve vector ix
% (identical equations and tolerances to one_run_full.m, single candidate)
P=D.P(:); B=D.B(:); IfP=D.IfP(:); IfB=D.IfB(:); nc=D.ncons; nR=D.nR;
BIDX=D.BIDX; maxd=size(BIDX,2); PS=PS(:); ix=ix(:);
Am=0.14*(ix==1)+13.5*(ix==2)+80*(ix==3); Bm=0.02*(ix==1)+1*(ix==2)+2*(ix==3);
PSMp=IfP(P)./(D.CT*PS(P)); PSMb=IfB./(D.CT*PS(B));
cp=Am(P)./(PSMp.^Bm(P)-1); cb=Am(B)./(PSMb.^Bm(B)-1);
vd=(PSMp>1.05)&(PSMb>1.05)&isfinite(cp)&isfinite(cb)&cp>0&cb>0;
cpv=cp; cpv(~vd)=0; cbv=cb; cbv(~vd)=1;
TM=D.TMSlb*ones(nR,1); req=-inf(nc+1,1);
for it=1:300
  rr=(D.CTI+cpv.*TM(P))./cbv; rr(~vd)=-inf; req(1:nc)=rr;
  tb=max(reshape(req(BIDX(:)),nR,maxd),[],2);
  Tn=max(TM,min(D.TMSub,tb));
  if max(abs(Tn-TM))<1e-10, TM=Tn; break; end
  TM=Tn;
end
vio=vd&((D.CTI+cpv.*TM(P))>cbv.*TM(B)+1e-6);
NV=sum(vio)+sum(~vd); Top=sum(cpv.*TM(P));
end
