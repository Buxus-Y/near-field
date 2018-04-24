%Integration based on the simpson rule and convergence criteria

function [Int_val] = specinte(w, n, d, min, max, ep_1, ep_2, err)
%used to calculate heat flux for each layer

         n1 = n;                                               
         [s_p,dkx] = func_p(w, n1, d, min, max, ep_1, ep_2);        
         temp1=(s_p(1)+4*sum(s_p(2:2:n1))+2*sum(s_p(3:2:n1-1))+s_p(n1+1))*dkx/3;

            if(temp1<=1e-50)
                temp2 = 0;
            else
                n1 = n1*2;
                [s_p,dkx] = func_p(w, n1, d, min, max, ep_1, ep_2);        
                temp2=(s_p(1)+4*sum(s_p(2:2:n1))+2*sum(s_p(3:2:n1-1))+s_p(n1+1))*dkx/3;
                while ((abs(temp2-temp1)/temp2) >= err)
                    temp1 = temp2;
                    n1 = n1*2;
                    [s_p,dkx] = func_p(w, n1, d, min, max, ep_1, ep_2);        
                    temp2=(s_p(1)+4*sum(s_p(2:2:n1))+2*sum(s_p(3:2:n1-1))+s_p(n1+1))*dkx/3;
                end
            end
                Int_val = temp2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s funtion divided by pi^2 in beta(kz) space for propagation waves
% 1--source; 2--receiver; 3--vacuum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s_p,dkx] = func_p(w, n, d, min, max,ep,ls)
 %spectral and wavevector dependent heat flux, propogating modes           
            
global c0
dkx = (max-min)/n;                                         
kx = zeros(n+1);        %horizontal component of wavevector k
s_p = zeros(n+1);
kz1=zeros(n+1);
kz2=zeros(n+1);
for ind=1:n+1
    kx(ind)= min+(ind-1)*dkx;
    kzs(ind)=sqrt(ep(ls)*w*w/(c0*c0)-kx(ind)^2);
for epind=1:length(ep)-1
  %  vertical component of wavevector k
      kz1(ind)=sqrt(ep(epind)*w*w/(c0*c0)-kx(ind)^2);
      kz2(ind)=sqrt(ep(epind+1)*w*w/(c0*c0)-kx(ind)^2);
       %   caculate frensnel coefficients of each layer
       rs=zeros(length(ep)-1);
       ts=zeros(length(ep)-1);
[ rs(epind),ts(epind)] = f_coefficients(ep(epind),ep(epind+1),kz1(ind),kz2(ind));
%caculate S-matrix
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
    end
    if layer>ls-1%for layer between soure and N
        Ss(layer+1,1,1)=(Ss(layer,1,1)*ts(layer)*exp(1i*kz1(ind)*d(layer)))/(1-Ss(layer,1,2)*rs(layer)*exp(2*1i*kz1(ind)*d(layer)));
        Ss(layer+1,1,2)=(Ss(layer,1,2)*exp(2*1i*kz1(ind)*d(layer)-ts(layer)))/(1-Ss(layer,1,2)*rs(layer)*exp(2*1i*kz1(ind)*d(layer)));
        Ss(layer+1,2,1)=(Ss(layer,1,1)*Ss(layer,2,2)*rs(layer)*exp(1i*kz1(ind)*d(layer)))/ts(layer)+Ss(layer,2,1);
        Ss(layer+1,2,2)=(Ss(layer,2,2)*(rs(layer)*Ss(layer,1,2)+1)*exp(1i*kz1(ind)*d(layer)))/ts(layer);
    end
end
    %caculate A B C and D
   coef=zeros(length(ep),4); %initialize A B C D matrix
   coef(ls,
   
end
    
