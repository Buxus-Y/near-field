%Integration based on the simpson rule and convergence criteria

function [Int_val] = dygreen_int(w, n, d, min, max,ep,ls,lc,err)
%used to calculate heat flux for each layer

         n1 = n;                                               
         [dgf,dkx] = dygreen(w, n, d, min, max,ep,ls,lc,err);        
         temp1=(dgf(1)+4*sum(dgf(2:2:n1))+2*sum(dgf(3:2:n1-1))+dgf(n1+1))*dkx/3;

            if(temp1<=1e-50)
                temp2 = 0;
            else
                n1 = n1*2;
                [dgf,dkx] = dygreen(w, n, d, min, max,ep,ls,lc,err);        
                temp2=(dgf(1)+4*sum(dgf(2:2:n1))+2*sum(dgf(3:2:n1-1))+dgf(n1+1))*dkx/3;
                while ((abs(temp2-temp1)/temp2) >= err)
                    temp1 = temp2;
                    n1 = n1*2;
                    [dgf,dkx] = dygreen(w, n, d, min, max,ep,ls,lc,err);        
                    temp2=(dgf(1)+4*sum(dgf(2:2:n1))+2*sum(dgf(3:2:n1-1))+dgf(n1+1))*dkx/3;
                end
            end
                Int_val = temp2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dgf,dkx] = dygreen(w, n, d, min, max,ep,ls,lc,errz)
 %spectral and wavevector dependent heat flux, TE mode           
            11;
global c0
dkx = (max-min)/n;                                         
kx = zeros(n+1);        %horizontal component of wavevector k
dgf = zeros(n+1);
kz1=zeros(n+1);
kz2=zeros(n+1);
for ind=1:n+1
    zn=n;
    dzx=d(ls)/zn;
for zx=sum(d(1:d(ls-1))):dzx:sum(d(1:ls))
    kx(ind)= min+(ind-1)*dkx;
	ks=ep(ls)*w*w/(c0*c0);
	kl=ep(lc)*w*w/(c0*c0);
    kzs(ind)=sqrt(ep(ls)*w*w/(c0*c0)-kx(ind)^2);
	kzl(ind)=sqrt(ep(lc)*w*w/(c0*c0)-kx(ind)^2);
	lep=length(ep);
for epind=1:lep-1
  %  vertical component of wavevector k
      kz1(ind)=sqrt(ep(epind)*w*w/(c0*c0)-kx(ind)^2);
      kz2(ind)=sqrt(ep(epind+1)*w*w/(c0*c0)-kx(ind)^2);
       %   caculate frensnel coefficients of each layer
       rs=zeros(length(ep)-1);
       ts=zeros(length(ep)-1);
[ rs(epind),ts(epind),rp(epind),tp(epind),] = f_coefficients(ep(epind),ep(epind+1),kz1(ind),kz2(ind));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calclulate for TE mode
%calclulate S-matrix
    S0=zeros(7,2,2);  % initialize S-matrix
    S0(1,:,:)=eye(2);
    Ss=S0;
    Ss(1,:,:)=eye(2);
    layer=epind;
    if layer<=ls-1 %for layer between 0 and soure
        S0(layer+1,1,1)=(S0(layer,1,1)*ts(layer)*exp(1i*kz1(ind)*d(layer)))/(1-S0(layer,1,2)*rs(layer)*exp(2*1i*kz1(ind)*d(layer)));
        S0(layer+1,1,2)=(S0(layer,1,2)*exp(2*1i*kz1(ind)*d(layer)-ts(layer)))/(1-S0(layer,1,2)*rs(layer)*exp(2*1i*kz1(ind)*d(layer)));
        S0(layer+1,2,1)=(S0(layer,1,1)*S0(layer,2,2)*rs(layer)*exp(1i*kz1(ind)*d(layer)))/ts(layer)+S0(layer,2,1);
        S0(layer+1,2,2)=(S0(layer,2,2)*(rs(layer)*S0(layer,1,2)+1)*exp(1i*kz1(ind)*d(layer)))/ts(layer);
  
    else %for layer between soure and N
        Ss(layer+1,1,1)=(Ss(layer,1,1)*ts(layer)*exp(1i*kz1(ind)*d(layer)))/(1-Ss(layer,1,2)*rs(layer)*exp(2*1i*kz1(ind)*d(layer)));
        Ss(layer+1,1,2)=(Ss(layer,1,2)*exp(2*1i*kz1(ind)*d(layer)-ts(layer)))/(1-Ss(layer,1,2)*rs(layer)*exp(2*1i*kz1(ind)*d(layer)));
        Ss(layer+1,2,1)=(Ss(layer,1,1)*Ss(layer,2,2)*rs(layer)*exp(1i*kz1(ind)*d(layer)))/ts(layer)+Ss(layer,2,1);
        Ss(layer+1,2,2)=(Ss(layer,2,2)*(rs(layer)*Ss(layer,1,2)+1)*exp(1i*kz1(ind)*d(layer)))/ts(layer);
    end
end
    %calculate A B C and D 
	ATE=zeros(lep);
	BTE=zeros(lep);
	CTE=zeros(lep);
	DTE=zeros(lep); %initialize A B C D matrix    coef=zeros(length(ep),4);
   Sp=exp(1i*kzs(ind)*(sum(d(1:ls-1)))-zx);% account for source S+
   %calculate  AL BL
   BTE(ls)=Ss(lep,2,1)*Sp/(1-Ss(lep,2,1)*S0(ls,1,2)); %caculate As B0 AN BS
   BTE(1)=S0(ls,2,2)*BTE(ls);
   ATE(ls)=S0(ls,1,2)*BTE(ls);
   ATE(lep)=Ss(lep,1,1)*(ATE(ls)+Sp);
   if lc<ls
   BTE(lc)=BTE(1)/S0(lc,1,1);
   ATE(lc)=S0(lc,1,2)*BTE(lc);
   else
   BTE(lc)=BTE(ls)-Ss(lc,2,1)*(ATE(ls)+Sp)/Ss(lc,2,2);
   ATE(lc)=Ss(lc,1,1)*(ATE(ls)+Sp)+Ss(lc,1,2)*BTE(lc);
   end
   Sd=exp(1i*kzs(ind)*(-sum(d(1:ls-1))));% account for source S-
    %calculate  CL DL
	CTE(ls)=S0(ls,1,2)*Sd/(1-S0(ls,1,2)*Ss(lep,2,1));
	DTE(ls)=Ss(lep,2,1)*CTE(ls);
	DTE(1)=S0(ls,2,2)*(DTE(ls)+Sd);
	CTE(lep)=Ss(lep,1,1)*CTE(ls);
   if lc<ls
   DTE(lc)=DTE(1)/S0(lc,2,2);
   CTE(lc)=S0(lc,1,2)*DTE(lc);
   else
   DTE(lc)=(BTE(ls)-Ss(lc,2,1)*CTE(ls))/Ss(lc,2,2);
   CTE(lc)=Ss(lc,1,1)*CTE(ls)+Ss(lc,1,2)*DTE(lc);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %calculate  for TMmode
     Sp0=zeros(7,2,2);  % initialize S-matrix
    Sp0(1,:,:)=eye(2);
    Sps=Sp0;
    Sps(1,:,:)=eye(2);
    layer=epind;
    if layer<=ls-1 %for layer between 0 and soure
        Sp0(layer+1,1,1)=(Sp0(layer,1,1)*tp(layer)*exp(1i*kz1(ind)*d(layer)))/(1-Sp0(layer,1,2)*rp(layer)*exp(2*1i*kz1(ind)*d(layer)));
        Sp0(layer+1,1,2)=(Sp0(layer,1,2)*exp(2*1i*kz1(ind)*d(layer)-tp(layer)))/(1-Sp0(layer,1,2)*rp(layer)*exp(2*1i*kz1(ind)*d(layer)));
        Sp0(layer+1,2,1)=(Sp0(layer,1,1)*Sp0(layer,2,2)*rp(layer)*exp(1i*kz1(ind)*d(layer)))/tp(layer)+Sp0(layer,2,1);
        Sp0(layer+1,2,2)=(Sp0(layer,2,2)*(rp(layer)*Sp0(layer,1,2)+1)*exp(1i*kz1(ind)*d(layer)))/tp(layer);
  
    else %for layer between soure and N
        Sps(layer+1,1,1)=(Sps(layer,1,1)*tp(layer)*exp(1i*kz1(ind)*d(layer)))/(1-Sps(layer,1,2)*rp(layer)*exp(2*1i*kz1(ind)*d(layer)));
        Sps(layer+1,1,2)=(Sps(layer,1,2)*exp(2*1i*kz1(ind)*d(layer)-tp(layer)))/(1-Sps(layer,1,2)*rp(layer)*exp(2*1i*kz1(ind)*d(layer)));
        Sps(layer+1,2,1)=(Sps(layer,1,1)*Sps(layer,2,2)*rp(layer)*exp(1i*kz1(ind)*d(layer)))/tp(layer)+Sps(layer,2,1);
        Sps(layer+1,2,2)=(Sps(layer,2,2)*(rp(layer)*Sps(layer,1,2)+1)*exp(1i*kz1(ind)*d(layer)))/tp(layer);
    end
    %calculATM A B C and D for TEmode
	ATM=zeros(lep);
	BTM=zeros(lep);
	CTM=zeros(lep);
	DTM=zeros(lep); %initialize A B C D matrix    coef=zeros(length(ep),4);
   Sp=exp(1i*kzs(ind)*(sum(d(1:ls-1))));% account for source S+
   %calculATM  AL BL
   BTM(ls)=Sps(lep,2,1)*Sp/(1-Sps(lep,2,1)*Sp0(ls,1,2)); %caculATM As B0 AN BS
   BTM(1)=Sp0(ls,2,2)*BTM(ls);
   ATM(ls)=Sp0(ls,1,2)*BTM(ls);
   ATM(lep)=Sps(lep,1,1)*(ATM(ls)+Sp);
   if lc<ls
   BTM(lc)=BTM(1)/Sp0(lc,1,1);
   ATM(lc)=Sp0(lc,1,2)*BTM(lc);
   else
   BTM(lc)=(BTM(ls)-Sps(lc,2,1)*(ATM(ls)+Sp))/Sps(lc,2,2);
   ATM(lc)=Sps(lc,1,1)*(ATM(ls)+Sp)+Sps(lc,1,2)*BTM(lc);
   end
   Sd=exp(1i*kzs(ind)*(zx-sum(d(1:ls-1))));% account for source S-
    %calculATM  CL DL
	CTM(ls)=Sp0(ls,1,2)*Sd/(1-Sp0(ls,1,2)*Sps(lep,2,1));
	DTM(ls)=Sps(lep,2,1)*CTM(ls);
	DTM(1)=Sp0(ls,2,2)*(DTM(ls)+Sd);
	CTM(lep)=Sps(lep,1,1)*CTM(ls);
   if lc<ls
   DTM(lc)=DTM(1)/Sp0(lc,2,2);
   CTM(lc)=Sp0(lc,1,2)*DTM(lc);
   else
   DTM(lc)=(BTM(ls)-Sps(lc,2,1)*CTM(ls))/Sps(lc,2,2);
   CTM(lc)=Sps(lc,1,1)*CTM(ls)+Sps(lc,1,2)*DTM(lc);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %calculate DGF
   gslxxE=1i*kzl(ind)/(2*ks*kl)*((ATM(lc)-BTM(lc))*exp(1i*-kzs(ind)*zx)+(-CTM(ls)+DTM(ls))*exp(1i*kzs(ind)*zx));
   gslyxH=kl/(2*ks)*((-ATM(lc)-BTM(lc))*exp(1i*-kzs(ind)*zx)+(CTM(lc)+DTM(lc))*exp(1i*kzs(ind)*zx));
   gslxzE=1i*kzl(ind)*kx(ind)/(2*kzs(ind)*ks*kl)*((-ATM(lc)+BTM(lc))*exp(1i*-kzs(ind)*zx)+(-CTM(ls)+DTM(ls))*exp(1i*kzs(ind)*zx));
   gslyzH=kl*kx(ind)/(2*ks*kzs(ind))*((ATM(lc)+BTM(lc))*exp(1i*-kzs(ind)*zx)+(CTM(ls)+DTM(ls))*exp(1i*kzs(ind)*zx));
   gslyyE=1i/(2*kzs(ind))*((ATE(lc)+BTE(lc))*exp(1i*-kzs(ind)*zx)+(CTE(ls)+DTE(ls))*exp(1i*kzs(ind)*zx));
   gslxyH=kzl(ind)/(2*kzs(ind))*((ATE(lc)-BTE(lc))*exp(1i*-kzs(ind)*zx)+(CTE(ls)-DTE(ls))*exp(1i*kzs(ind)*zx));
   dgfz=gslxxE*conj(gslyxH)+gslxzE*conj(gslyzH)+gslyyE*conj(gslxyH);
end
text1=(dgfz(1)+4*sum(dgfz(2:2:zn))+2*sum(dgfz(3:2:zn-1))+dgfz(zn+1))*dzx/3;
   if(text1<=1e-50)
                text2= 0;
            else
                zn = zn*2;       
                text2=(dgfz(1)+4*sum(dgfz(2:2:zn))+2*sum(dgfz(3:2:zn-1))+dgfz(zn+1))*dzx/3;
                while ((abs(text2-text1)/text2) >= errz)
                    text1 = text2;
                    zn = zn*2;       
                    text2=(dgfz(1)+4*sum(dgfz(2:2:zn))+2*sum(dgfz(3:2:zn-1))+dgfz(zn+1))*dzx/3;
                end
   end
   dgf=text2;
end
end
end

    
