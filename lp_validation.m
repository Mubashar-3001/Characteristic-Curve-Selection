% lp_validation.m -- inner Bellman-Ford relaxation vs linprog, 62-pair set, uniform SI
D=load('docr_corrected_data62.mat');
nR=D.nR; P=D.P(:); B=D.B(:); IfP=D.IfP(:); IfB=D.IfB(:); CT=D.CT; CTI=D.CTI;
nc=D.ncons; BIDX=D.BIDX; TMSlb=D.TMSlb; TMSub=D.TMSub; Ac=0.14; Bc=0.02; PSMmin=1.05;
ubr=D.PSub*ones(nR,1);
for r=1:nR, lim=D.PSub; if any(P==r)&&IfP(r)>0, lim=min(lim,IfP(r)/(CT*PSMmin)); end
  kk=find(B==r); for k=kk(:)', lim=min(lim,IfB(k)/(CT*PSMmin)); end; ubr(r)=max(lim,D.PSlb); end
rng(4242,'twister'); nTest=40; res=zeros(nTest,6);
for t=1:nTest
  PS=D.PSlb+(ubr-D.PSlb).*rand(nR,1);
  gp=Ac./((IfP(P)./(CT*PS(P))).^Bc-1); gb=Ac./((IfB./(CT*PS(B))).^Bc-1);
  TMS=TMSlb*ones(nR,1);
  for s=1:300, old=TMS; rq=-inf(nc+1,1); rq(1:nc)=(CTI+gp.*TMS(P))./gb;
    TMS=max(TMS,min(TMSub,max(rq(BIDX),[],2))); if all(abs(TMS-old)<1e-12), break; end; end
  f=zeros(nR,1); for k=1:nc, f(P(k))=f(P(k))+gp(k); end
  A=zeros(nc,nR); b=-CTI*ones(nc,1);
  for k=1:nc, A(k,P(k))=A(k,P(k))+gp(k); A(k,B(k))=A(k,B(k))-gb(k); end
  [~,fv,flag]=linprog(f,A,b,[],[],TMSlb*ones(nR,1),TMSub*ones(nR,1),optimoptions('linprog','Display','none'));
  if isempty(fv), fv=NaN; end
  res(t,:)=[t sum(gp.*TMS(P)) fv flag s max(A*TMS-b)];
end
ok=res(:,4)==1; rd=100*(res(ok,2)-res(ok,3))./res(ok,3);
LPV=struct('res',res,'nFeas',sum(ok),'nInf',sum(~ok),'maxRel',max(abs(rd)),'medRel',median(abs(rd)), ...
  'maxResFeas',max(res(ok,6)),'maxResInf',max(res(~ok,6)),'swMin',min(res(:,5)),'swMax',max(res(:,5)),'ubr',ubr);
save('lp_validation62.mat','LPV');