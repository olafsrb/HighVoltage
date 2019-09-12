clear all;
close all;
format long;
clc;
%---------------------------------------------------------------------------------------------------%
                                           %CONSTANTS%
%---------------------------------------------------------------------------------------------------%

n_iter = 10000;
aa = 1.5;
k= 1;
A= 1.5; % Ratio for H and W later
%Nr of grids
nx=200*k; %lines
ny=140*k;%colums

CenterLineRow = nx/2;
CenterLineCol = ny/2;

RatioBoundNX = 8;
RatioBoundNY = 8;

RowOffCenter = 0.5*(nx/RatioBoundNX);
ColOffCenter = 0.5*(ny/RatioBoundNY);

nxa = round(CenterLineRow-RowOffCenter);
nxb = round(CenterLineRow + RowOffCenter);
nya = round(CenterLineCol-ColOffCenter);
nyb = round(CenterLineCol+ColOffCenter);

%put k constant to follow the grid, these are like measurment probes
iprobe_x= 5*k;
iprobe_y= 5*k;
%permativity constant for epoxy
Eps_r1 = 10;  
%permativity const for air
Eps_r2 = 1;
%permativity Boundaries
Eps_r3 = (Eps_r1+Eps_r2)/2;
%---------------------------------------------------------------------------------------------------%
                                        %Voltage,eps and Marker Matrices%
%---------------------------------------------------------------------------------------------------%

% voltage potentional matrix
u=0.5*ones(nx,ny);
%u(nxa:nxb,nya:nyb)=1;
 % u(:,1) = 0;
 % u(:,ny) = 0;
 u(1,:) = 1;
 u(nx,:)= 0;
 
%Markers
iu=zeros(nx,ny);

iu(nxa,nya:nyb) = -1;

iu(nxb,nya:nyb) = -1;

iu(nxa:nxb,nya) = -1;

iu(nxa:nxb,nyb) = -1;

%Make a eps matrix
Eps = ones(nx,ny)*Eps_r1; %epoxy
Eps(nxa:nxb,nya:nyb)= Eps_r3; %air+epoxy
Eps(nxa+1:nxb-1,nya+1:nyb-1) = Eps_r2; %air
%---------------------------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------------------------%
t=cputime; 
it = 1;
tol = 1;
u0 = 1;
while (tol >= 1e-6)  
    for ix = 2:nx-1
        for iy = 2:ny-1
           if ( iu(ix,iy) == 0) 
                 rnm = u(ix,iy)- ( u(ix+1,iy)+ u(ix-1,iy) + u(ix,iy+1) + u(ix,iy-1) )/4;
                 
                 u(ix,iy) = u(ix,iy) - aa*rnm;
                 
           end 
            if ( iu(ix,iy) == -1)
                u(ix,iy) = ( (Eps(ix+1,iy))*(u(ix+1,iy)) + (Eps(ix-1,iy))*(u(ix-1,iy)) + (Eps(ix,iy+1))*(u(ix,iy+1)) + (Eps(ix,iy-1))*(u(ix,iy-1)) )/( Eps(ix+1,iy)+ Eps(ix-1,iy)+ Eps(ix,iy+1)+ Eps(ix,iy-1)); %needs to be corrected to our problem at the boundaries
                   
            end
  
        end
    end
      u(2:nx-1,1) = u(2:nx-1,2);
      u(2:nx-1,ny) = u(2:nx-1,ny-1); 
   tol = abs(u(iprobe_x,iprobe_y)-u0);
   u0 = u(iprobe_x,iprobe_y);
   res_n(it,1)= tol;
   it = it+1;
end;

cputime-t;
save potential u iu
%---------------------------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------------------------%
% 
figure(1)
spy(iu)
grid on
%set(gca,'xlim',[0,ny],'ylim',[0,nx],'ydir','reverse','GridLineStyle','none','plotboxaspectratio',[nx+1 ny+1 1]);
set(gca, 'GridLineStyle', ':')
 figure(2)
 plot(log10(abs(res_n)),'r.');
 grid on;
