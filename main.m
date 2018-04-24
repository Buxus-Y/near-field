clear all
clc
tic

% Physical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 1.054571596e-34;           %Planck's constant over 2*pi
c0 = 2.99792458e+8;            %speed of light in vacuum
kb = 1.3806503e-23;
qe = 1.602176462e-19;          %elementary electric charge %%for doped silicon calculation
e0 = 8.854187817e-12;          %electric permittivity %% for doped silicon calculation
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
T1 = 300;                      %temperature of the emitter  @@@@@@@@
T2 = 290;                      %temperature of the receiver@@@@@@@@
d = 1e-8;%gap spacing
lp=1e-7;
ln=dp;
li=1e-6;
lr=5e-5;
darray=[0,lp,li,ln,d,lr];
ep=[epA,epn,epi,epp,1,epr,epA];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiaization of different parameters
w1=2.0e12; %rad/s about 940 um in wavelength
w2=2.0e15; %rad/s about 0.94 um in wavelength
dw=2.0e12;
nw=floor((w2-w1)/dw+1)  %1001 2001
             a0 = 0;
             nk = 1001; %It should be sufficiently large to obtain good accuracy
             const1 = 3;
	     const2 = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculation of total energy flux
      for ii=1:nw
        for j=1
             w(ii)=w1+(ii-1)*dw;
             %dielectric function
             epsi1=Lorentz_GaAs(w(ii));%for emitter
             epsi2=Drude_Si(w(ii));%for receiver
% Calculating the radiative heat transfer
%------------------------------------------
             a1 = (w(ii)-10)/c0;    %a0 and a1 for the integration of propagating waves. 
             ae1 = (w(ii)+10)/c0;     %integration limit of evanescient waves
             ae2 = max(0.6/d,const1*w(ii)/c0); 
             ae3 = max(5/d,const2*w(ii)/c0);
             
             Q1(ii)=h*w(ii)/(exp(h*w(ii)/(kb*T1))-1);    %% Planck Oscillator
             Q2(ii)=h*w(ii)/(exp(h*w(ii)/(kb*T2))-1);
             [T_w_p] = Uni_spec(w(ii), nk, d, a0, a1, epsi_O1, epsi_E1, epsi_O2, epsi_E2);
             [T_w_e1]= Uni_spec(w(ii), nk, d, ae1, ae2, epsi_O1, epsi_E1, epsi_O2, epsi_E2);
             [T_w_e2]= Uni_spec(w(ii), nk, d, ae2, ae3, epsi_O1, epsi_E1, epsi_O2, epsi_E2);
             S_weach(ii,j) = (T_w_p+T_w_e1+T_w_e2)*(Q1(ii)-Q2(ii));

%Note: S_w gives an array of the spectral heat flux@@@@@@@@@@@@@
    
        end
