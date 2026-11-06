function [gT,gN,gTM,gIX,conv,LOG] = one_run_full_logged(seed,N,T,ri,LB,UB,rg,a1,a2,a3,a4,nR,P,B,IfP,IfB,CT,CTI,nc,BIDX,maxd,TMSlb,TMSub,EI)
vmax=0.5*rg;
rng(seed,'twister');
X=LB+rg*rand(N,nR); V=zeros(N,nR); V1=V;V2=V;V3=V;
gP=inf; gX=X(1,:); gT=inf; gN=inf; gTM=[]; gIX=[]; conv=nan(1,T+1);
LOG.H=nan(1,T); LOG.reinj=false(1,T); LOG.sweeps=zeros(1,T+1); LOG.hitcap=false(1,T+1); LOG.nviol_pop=zeros(1,T+1);
for t=0:T
  PSt=X'; PSMr=IfP./(CT*PSt); bad=PSMr<=1.05;
  switch ri
    case 1, ix=ones(nR,N);
    case 2, ix=2*ones(nR,N);
    case 3, ix=3*ones(nR,N);
    case 4, ix=ones(nR,N); ix(EI,:)=3;
    case 5, ix=ones(nR,N); ix(PSMr>5&PSMr<=20)=2; ix(PSMr>20)=3;
    case 6
      cS=0.14./(PSMr.^0.02-1); cV=13.5./(PSMr-1); cE=80./(PSMr.^2-1);
      cS(bad)=inf; cV(bad)=inf; cE(bad)=inf;
      [~,ix]=min(cat(3,cS,cV,cE),[],3);
  end
  Am=0.14*(ix==1)+13.5*(ix==2)+80*(ix==3);
  Bm=0.02*(ix==1)+1*(ix==2)+2*(ix==3);
  PSMp=IfP(P)./(CT*PSt(P,:)); PSMb=IfB./(CT*PSt(B,:));
  cp=Am(P,:)./(PSMp.^Bm(P,:)-1); cb=Am(B,:)./(PSMb.^Bm(B,:)-1);
  vd=(PSMp>1.05)&(PSMb>1.05)&isfinite(cp)&isfinite(cb)&cp>0&cb>0;
  cpv=cp; cpv(~vd)=0; cbv=cb; cbv(~vd)=1;
  TM=TMSlb*ones(nR,N); req=-inf(nc+1,N);
  conv_ok=false;
  for it=1:300
    rr=(CTI+cpv.*TM(P,:))./cbv; rr(~vd)=-inf; req(1:nc,:)=rr;
    tb=reshape(max(reshape(req(BIDX(:),:),nR,maxd,N),[],2),nR,N);
    Tn=max(TM,min(TMSub,tb));
    if max(max(abs(Tn-TM)))<1e-10, TM=Tn; conv_ok=true; break; end
    TM=Tn;
  end
  vio=vd&((CTI+cpv.*TM(P,:))>cbv.*TM(B,:)+1e-6);
  NV=sum(vio,1)'+sum(~vd,1)'; Topv=sum(cpv.*TM(P,:),1)';
  Phi=NV*1e6+Topv;
  LOG.sweeps(t+1)=it; LOG.hitcap(t+1)=~conv_ok; LOG.nviol_pop(t+1)=sum(NV>0);
  [mP,mi]=min(Phi);
  if mP<gP, gP=mP; gX=X(mi,:); gT=Topv(mi); gN=NV(mi); gTM=TM(:,mi); gIX=ix(:,mi); end
  conv(t+1)=gT;
  if t==T, break; end
  fw=max(Phi); fb=min(Phi); mM=(fw-Phi)/(fw-fb+eps); Mm=mM/(sum(mM)+eps);
  G=100*exp(-20*t/T); Kb=max(1,round(N-(N-1)*t/T));
  [~,ord]=sort(Phi,'ascend'); acc=zeros(N,nR);
  for jj=ord(1:Kb)'
    dX=X(jj,:)-X; Rr=sqrt(sum(dX.^2,2))+eps;
    acc=acc+rand(N,1).*(G*Mm(jj)./Rr).*dX;
  end
  d=sqrt(sum((X-gX).^2,2));
  if sum(d)<eps, Hn=0; else, pp=d/sum(d); Hn=-sum(pp.*log(pp+eps))/log(N); end
  LOG.H(t+1)=Hn;
  Vn=a1*V+a2*V1+a3*V2+a4*V3+2*rand(N,nR).*(1+Hn).*acc+2*rand(N,nR).*(gX-X);
  Vn=max(min(Vn,vmax),-vmax); V3=V2;V2=V1;V1=V;V=Vn;
  X=max(min(X+V,UB),LB);
  if Hn<0.5
    LOG.reinj(t+1)=true;
    nre=max(1,round(0.2*N)); [~,wo]=sort(Phi,'descend'); id=wo(1:nre);
    X(id,:)=LB+rg*rand(nre,nR); V(id,:)=0;V1(id,:)=0;V2(id,:)=0;V3(id,:)=0;
  end
end
end
