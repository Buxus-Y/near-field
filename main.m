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
bias=[0:0.001:1.3];
epsi=[epA,epn,epi,epp,1,epr,epA]; %dielectric function
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
             w(ii)=w1+(ii-1)*dw;
             
% Calculating the radiative heat transfer
%------------------------------------------
             a1 = (w(ii)-10)/c0;    %a0 and a1 for the integration of propagating waves. 
             ae1 = (w(ii)+10)/c0;     %integration limit of evanescient waves
             ae2 = max(0.6/d,const1*w(ii)/c0); 
             ae3 = max(5/d,const2*w(ii)/c0);
             errp=0.01;
			 err1=0.1;
			 err2=0.01;
             Q1(ii)=h*w(ii)/(exp(h*w(ii)/(kb*T1))-1);    %% Planck Oscillator
             Q2(ii)=h*w(ii)/(exp((h*w(ii)-qe*bias)/(kb*T1))-1);
			 Q3(ii)=h*w(ii)/(exp(h*w(ii)/(kb*T1))-1);
			 T_w=zeros(3,4);
             [T_w(1,1)] = dygreen_int(w(ii), nk, darray,a0, a1, epsi,2 , 6,errp);
			 [T_w(1,2)] = dygreen_int(w(ii), nk, darray,a0, a1, epsi,3 , 6,errp);
			 [T_w(1,3)] = dygreen_int(w(ii), nk, darray ,a0, a1, epsi,4 , 6,errp);
			 [T_w(1,4)] = dygreen_int(w(ii), nk, darray ,a0, a1, epsi,6 , 4,errp);
             [T_w(2,1)]= dygreen_int(w(ii), nk, darray, ae1, ae2, epsi, 2, 6, err1);
			 [T_w(2,2)]= dygreen_int(w(ii), nk, darray, ae1, ae2, epsi, 3, 6, err1);
			 [T_w(2,3)]= dygreen_int(w(ii), nk, darray, ae1, ae2, epsi, 4, 6, err1);
			 [T_w(2,4)]= dygreen_int(w(ii), nk, darray, ae1, ae2, epsi, 6, 4, err1);
             [T_w(3,1)]= dygreen_int(w(ii), nk, darray, ae2, ae3, epsi, 2, 6, err2);
			 [T_w(3,2)]= dygreen_int(w(ii), nk, darray, ae2, ae3, epsi, 3, 6, err2);
			 [T_w(3,3)]= dygreen_int(w(ii), nk, darray, ae2, ae3, epsi, 4, 6, err2);
			 [T_w(3,4)]= dygreen_int(w(ii), nk, darray, ae2, ae3, epsi, 6, 4, err2);
             S_w(ii) = sum(T_w(1,1:3))*Q1(ii)+sum(T_w(2,1:3))*Q2(ii)+sum(T_w(3,1:3))*Q1(ii)-sum(T_w(4,1:3))*Q3(ii);
%Note: S_w gives an array of the spectral heat flux@@@@@@@@@@@@@
    
      end
  
             Q_t=(S_w(1)+4*sum(S_w(2:2:(nw-1)))+2*sum(S_w(3:2:(nw-2)))+S_w(nw))*dw/3;