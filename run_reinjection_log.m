% run_reinjection_log.m -- instrumented campaign (identical RNG stream to one_run_full)
D=load('docr_corrected_data62.mat');
nR=D.nR; P=D.P; B=D.B; IfP=D.IfP(:); IfB=D.IfB(:); CT=D.CT; CTI=D.CTI;
nc=D.ncons; BIDX=D.BIDX; maxd=size(BIDX,2); TMSlb=D.TMSlb; TMSub=D.TMSub;
LB=D.PSlb; UB=D.PSub; rg=UB-LB; EI=D.EI;
al=0.1; a1=al; a2=al*(1-al)/2; a3=al*(1-al)*(2-al)/6; a4=al*(1-al)*(2-al)*(3-al)/24;
LOGS=struct();
for ri=1:6
  nre=zeros(50,1); swmean=nre; swmax=nre; caps=nre; tops=nre; Hend=nre;
  for k=1:50
    [gT,~,~,~,~,Lg]=one_run_full_logged(5000+k,50,300,ri,LB,UB,rg,a1,a2,a3,a4,nR,P,B,IfP,IfB,CT,CTI,nc,BIDX,maxd,TMSlb,TMSub,EI);
    tops(k)=gT; nre(k)=sum(Lg.reinj); swmean(k)=mean(Lg.sweeps); swmax(k)=max(Lg.sweeps); caps(k)=sum(Lg.hitcap); Hend(k)=Lg.H(end);
  end
  LOGS.(sprintf('R%d',ri))=struct('tops',tops,'nre',nre,'swmean',swmean,'swmax',swmax,'caps',caps,'Hend',Hend);
end
save('reinjection_log62.mat','LOGS');