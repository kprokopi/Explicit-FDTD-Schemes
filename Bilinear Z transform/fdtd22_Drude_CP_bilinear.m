%1D FDTD simulation
%A wave impinges upon a gold layer (20 nm thick)
%The gold is modelled by the Drude - critical point (DCP) model
%and the FDTD method is based on the bilinear Z-transform

% Prokopidis Kostas, 25/12/2012

clear all;
NZ=8000; %computational domain
m0=12.566E-7; %permability of free space
e0=8.854187817620E-12; %permittivity of free space
c0=1/sqrt(e0*m0); %velocity of light

%-----------------------------------------
Q=0.3; %Courant number
dz=1e-9; %space step
dt=Q*dz/c0; %time step

Nt=30000; %Number of time steps
%-----------------------------------------

%slab of 20nm
za=4000;
zb=4019;

% Parameters of Au from paper of Vial 2011
 epsilon_inf=1.1431;
 omega_D=1.3202e16;
 gamma=1.0805e14;
 
 A=[0.26698 3.0834];
 phi_m=[-1.2371 -1.0968];
 Omega=[3.8711e15 4.1684e15];
 Gamma=[4.4642e14 2.3555e15];

% Parameters of Au from paper of Vial 2007
% epsilon_inf=1.0300;
% omega_D=1.3064e16;
% gamma=1.1274e14;
% 
% A=[0.86822 1.37];
% phi_m=[-0.60756 -0.087341];
% Omega=[4.0812e15 6.4269e15];
% Gamma=[7.3277e14 6.7371e14];

%matrices
EX=zeros(3,NZ); 
HY=zeros(1,NZ-1);
DX=zeros(1,NZ);

%auxiliary variable
M_er=length(A)+1;
F=zeros(M_er,2,NZ-1);


% source constants
zsource=3; %the source location

%source (A)
%t0=100;
%func(1:Nt)=exp(-((1:Nt)*dt-3*t0*dt).^2./(t0*dt)^2); 

%source (B)
f_high=3e15;
tw=2/(pi*f_high);
t0=4*tw;
func(1:Nt)=exp(-((1:Nt)*dt-t0).^2./(tw)^2); 

%FDTD coefficients
a_0=2.*A.*Omega.*(Omega.*cos(phi_m)-Gamma.*sin(phi_m));
a_1=-2.*A.*Omega.*sin(phi_m);
b_0=Omega.^2+Gamma.^2;
b_1=2*Gamma;

%The Drude term is included here
a_0=[a_0 omega_D^2];
a_1=[a_1 0];
b_0=[b_0 0];
b_1=[b_1 gamma];

%Coefficients of the term in Z-domain
c_0=2.*a_1.*dt+a_0.*(dt)^2;
c_1=2.*a_0.*(dt)^2;
c_2=a_0.*(dt)^2-2.*a_1.*dt;

d_0=4+2.*b_1.*dt+b_0.*(dt)^2;
d_1=-8+2.*b_0.*(dt)^2;
d_2=4-2.*b_1.*dt+b_0.*(dt)^2;

c_EX=e0*(epsilon_inf+sum(c_0./d_0));


cH=dt/(m0*dz);
cE=dt/(e0*dz);

c_Mur=(c0*dt-dz)/(c0*dt+dz);

for n=1:Nt
        
   %Mur's ABC at k=NZ (part a)
   EX(1,NZ)=EX(1,NZ-1)-c_Mur.*EX(1,NZ);
   
   %Mur's ABC at k=1 (part a)
   EX(1,1)=EX(1,2)-c_Mur.*EX(1,1);
   
   %free space
   EX(1,2:za-1)=EX(1,2:za-1)-cE*(HY(2:za-1)-HY(1:za-2));
   
    
  %interface k=za
  %more code to be added here later ...
  

  % slab of gold [za,zb]
     
  %calculate DX
  DX(za:zb)=DX(za:zb)-(dt/dz)*(HY(za:zb)-HY(za-1:zb-1));
    
  %update EX
  for k=za:zb
      
  %back-storing EX
  EX(3,k)=EX(2,k);
  EX(2,k)=EX(1,k);
      
  temp=0;
    for p=1:M_er
        temp=temp+EX(2,k)*(c_1(p)/d_0(p))+EX(3,k)*(c_2(p)/d_0(p))-F(p,1,k)*(d_1(p)/d_0(p))-F(p,2,k)*(d_2(p)/d_0(p));
    
    end
   EX(1,k)=(DX(k)-e0*temp)/c_EX; 
    
   for p=1:M_er
       %back-storing F
   Temp=F(p,2,k);  
   F(p,2,k)=F(p,1,k); 
    
  %update F
 
  F(p,1,k)=(c_0(p)/d_0(p))*EX(1,k)+(c_1(p)/d_0(p))*EX(2,k)+(c_2(p)/d_0(p))*EX(3,k)...
      -(d_1(p)/d_0(p))*F(p,2,k)-(d_2(p)/d_0(p))*Temp;    
  end
   
  end  
  
  
  %interface k=zb
  %more code to be added here later...
 
   %free space
   EX(1,zb+1:NZ-1)=EX(1,zb+1:NZ-1)-cE*(HY(zb+1:NZ-1)-HY(zb:NZ-2));
      
   %Mur's ABC at k=NZ (part b)
   EX(1,NZ)=EX(1,NZ)+c_Mur.*EX(1,NZ-1);
   
   %Mur's ABC at k=1 (part b)
   EX(1,1)=EX(1,1)+c_Mur.*EX(1,2);
   
   %source (soft)
   EX(1,zsource)=EX(1,zsource)+func(n); 
             
   %Observation point
   
   PCP_Z(n)=EX(1,3999); 
   TCP_Z(n)=EX(1,4020);
    
   %H -component   
   HY(1:NZ-1)=HY(1:NZ-1)-cH*(EX(1,2:NZ)-EX(1,1:NZ-1));
   
   %movie
   if rem(n,100)==0;
       s=int2str(n);
       n2=n/100;
       clf;
       plot(EX(1,:));
       title(n);
       axis([0 8000 -1.5 1.5]);
       if n==100;
	M=moviein(100);
       end;
       hold;
       clc;
       M(:,n2)=getframe;
   end;  
   
   
   
end
%t=dt:dt:Nt*dt;
%plot(t./1e-9,PCP,'g');
%xlabel('t (nsec)');







