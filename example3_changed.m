clear all;
close all;
format long;
clc;
%---------------------------------------------------------------------------------------------------%
                                           %CONSTANTS%
%---------------------------------------------------------------------------------------------------%

n_iter = 10000;
aa = 1.5;
k=1;
A = 2;
%Nr of grids
nx=100*k*A; %lines
ny=140*k;%colums


nxa = ceil(nx /2.5); %this is boundary for PosZero in a line
nxb = ceil(nx/1.67); % this is boundary for PosFinal in line
nya = ceil(ny/2.4); % this is boundary for PosZero in a row
nyb = ceil(ny/1.71); % this is boundary for PosFinal in a row

%put k constant to follow the grid, these are like measurment probes
iprobe_x=5*k;
iprobe_y=5*k;
%permativity constant for epoxy
Eps_r1 = 10;  
%permativity const for air
Eps_r2 = 1;
Eps_r3 = (Eps_r1+Eps_r2)/2;
%---------------------------------------------------------------------------------------------------%
                                        %MATRIX%
%---------------------------------------------------------------------------------------------------%

%Make a eps matrix
Eps = ones(nx,ny)*Eps_r1;
% voltage potentional matrix
u=0.5*ones(nx,ny);
% u(nxa:nxb,nya:nyb)=1;
%  u(:,1) = 0;
%  u(:,ny) = 0;
 u(1,:) = 0;
 u(nx,:)= 0;
%---------------------------------------------------------------------------------------------------%
                                          %Markers%
%---------------------------------------------------------------------------------------------------%
iu=zeros(nx,ny);
% fylli hann a� innan
iu(nxa:nxb,nya:nyb)=1;
iu(nxa,nya:nyb) = -1;

iu(nxb,nya:nyb) = -1;

iu(nxa:nxb,nya) = -1;

iu(nxa:nxb,nyb) = -1;

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

figure(1)
spy(iu)
figure(2)
plot(log10(abs(res_n)),'r.');
grid on;