clc
clear all 
close all
%% Constants
% semi-angle at half power
theta=70;
%Lambertian order of emission
m=-log10(2)/log10(cosd(theta));
%transmitted optical power by individual LED
P_LED=20;
% number of LED array nLED*nLED
nLED=60;
%Total transmitted power
P_total=nLED*nLED*P_LED;
%detector physical area of a PD
Adet=1e-4;
%reflection coefficient
rho=0.8;
%gain of an optical filter; ignore if no filter is used
Ts=1;
%refractive index of a lens at a PD; ignore if no lens is used
index=1.5;
%FOV of a receiver
FOV=70;
%gain of an optical concentrator; ignore if no lens is used
G_Con=(index.^2)/(sind(FOV).^2);
%Number of RIS elements per row
M = 5; %max value 200 (200*0.025 = 5m - width of the room)
%Number of rows of RIS elements
N = 5; %max value 86 (86*0.025 = 2.15m - height of the room)
%Carrier Wavelength
lambda = 0.1; %In meters
%The width and height of an RIS element
w = lambda/4; % 0.025m
%Separation distance between adjacent RIS elements
dsep = lambda/2;
%% Room layout with transmitter and receiver position
% room dimension in meter
lx=5; ly=5; lz=2.15; %standard ceiling height 2.75m
% position of Transmitter (LED);
[XT,YT,ZT]=meshgrid(0,0,lz);
% number of grid in each surface
Nx=lx*2; Ny=ly*2; Nz=round(lz*2);
% calculation grid area
dA=lz*ly/(Ny*N);
x=linspace(-lx/2,lx/2,Nx);
xRIS = linspace(-lx/2,lx/2,M);
y=linspace(-ly/2,ly/2,Ny);
z=linspace(-lz/2,lz/2,Nz);
zRIS = linspace(-lz/2,lz/2,N);
[XR,YR,ZR]=meshgrid(x,y,-lz/2);
%% Power Calculation
% %% Channel Gain due to reflection
SUM_Power = zeros(Nx,Ny);
for xx=1:length(XT)
    for yy = 1:length(YT)
        if lx>=M*w && lz>=N*w
            Separate_Source_Power= Power_CAL_RIS(xx,yy,m,P_total,Adet,rho,Ts,FOV,G_Con,lx,ly,lz,XT,YT,ZT,lambda,dsep,M,N,Nx,Ny,Nz,dA,xRIS,zRIS,x,y,z);
            %surfc(x,y,Separate_Soucer_Power);   % This commands to draw the
            %reflection of each source individually. 
            %  hold on   
            SUM_Power = SUM_Power + Separate_Source_Power;
        else 
            disp('RIS dimension exceeds room dimension')
            return
        end 
    end
end
P_rec_total_1ref_ris=SUM_Power;%P_rec_A1+P_rec_A2+P_rec_A3+P_rec_A4;
P_rec_1ref_dBm=10*log10(P_rec_total_1ref_ris);

%% Directed Path
[XTd,YTd]= meshgrid(0,0);
xd=linspace(-lx/2,lx/2,Nx); %This is direct path, so receiver on the ground; Nx won't be replaced with M
yd=linspace(-ly/2,ly/2,Ny);
[XRd,YRd]=meshgrid(xd,yd); %Receiver x-y plane (on the floor) 
SUM_DIR_POWER = zeros(Nx,Ny);
   inn=1;
for xxd = 1 :length(XTd)
    for yyd = 1 :length(YTd)
     
  Direct_Power (:,:,inn) =  Direct_SUM_POWER(xxd,yyd,m,P_total,Adet,Ts,FOV,G_Con,lz,XTd,YTd,XRd,YRd);
%   % The below commands draw the direct power of each source individually.  
%   surfc(x,y,10*log10(Direct_Power(:,:,inn)));
%   xlabel('X (m)');
%   ylabel('Y (m)');
%   zlabel('Received power (dBm)');
%   axis([-lx/2 lx/2 -ly/2 ly/2 min(min(10*log10(Direct_Power(:,:,inn)))) 2*max(max(10*log10(Direct_Power(:,:,inn))))]);
%   hold on
%   title('Directed Path Received Power for individual sources');
    SUM_DIR_POWER = SUM_DIR_POWER + Direct_Power(:,:,inn);
      inn= inn+1;
    end
end
  P_rec_totald= SUM_DIR_POWER;%P_rec_A1d+P_rec_A2d+P_rec_A3d+P_rec_A4d;
P_rec_dBmd=10*log10(P_rec_totald);
%% Plotting Directed
% figure
% surfc(x,y,P_rec_dBmd);
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('Received power (dBm)');
% axis([-lx/2 lx/2 -ly/2 ly/2 min(min(P_rec_dBmd)) max(max(P_rec_dBmd))]);
% title('Directed Path Received Power');
% figure 
% contour(x,y,P_rec_dBmd);
% hold on
% title('Directed Path Received Power');
%mesh(x,y,P_rec_dBm);
%% Plotting First Reflection
% figure
% surfc(x,y,abs(P_rec_total_1ref_ris));
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('Received power');
% axis([-lx/2 lx/2 -ly/2 ly/2 min(min(abs(P_rec_total_1ref_ris))) max(max(abs(P_rec_total_1ref_ris)))]);
% title('First Reflection Received Power');
% figure
% surfc(x,y,abs(P_rec_1ref_dBm));
% xlabel('X (m)');
% ylabel('Y (m)');
% zlabel('Received power (dBm)');
% axis([-lx/2 lx/2 -ly/2 ly/2 min(min(abs(P_rec_1ref_dBm))) max(max(abs(P_rec_1ref_dBm)))]);
% title('First Reflection Received Power in dBm');
%% Plotting Total
P_total=P_rec_total_1ref_ris+P_rec_totald;
P_total_dBm=10*log10(P_total);
figure
surfc(x,y,abs(P_total));
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Received power');
axis([-lx/2 lx/2 -ly/2 ly/2 min(min(abs(P_total))) max(max(abs(P_total)))]);
title('Total Received Power');
%% Plotting Total in dBm
figure
surfc(x,y,abs(P_total_dBm));
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Received power (dBm)');
axis([-lx/2 lx/2 -ly/2 ly/2 min(min(abs(P_total_dBm))) max(max(abs(P_total_dBm)))]);
title('Total Received Power in dBm');