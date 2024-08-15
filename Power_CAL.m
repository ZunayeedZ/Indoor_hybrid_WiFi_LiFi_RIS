function [Powers] = Power_CAL(xx,yy,m,P_total,Adet,rho,Ts,FOV,G_Con,lx,ly,lz,XT,YT,ZT,Nx,Ny,Nz,dA,x,y,z)
%first transmitter calculation
TP1=[XT(1,xx,1) YT(yy,1,1) ZT(1,1,1)];
% transmitter position
TPV=[0 0 -1];
% transmitter position vector
RPV=[0 0 1];
% receiver position vector
WPV1=[1 0 0];
% position vector for wall 1
for ii=1:Nx
for jj=1:Ny
RP=[x(ii) y(jj) -lz/2];
% receiver position vector
h1(ii,jj)=0;
% reflection from North face
for kk=1:Ny
for ll=1:Nz
WP1=[-lx/2 y(kk) z(ll)];
% point of incidence in wall
D1=sqrt(dot(TP1-WP1,TP1-WP1));
cos_phi= abs(WP1(3)- TP1(3))/D1;
cos_alpha=abs(TP1(1)- WP1(1))/D1;
D2=sqrt(dot(WP1-RP,WP1-RP));
cos_beta=abs(WP1(1)- RP(1))/D2;
cos_psi=abs(WP1(3)- RP(3))/D2;
if abs(acosd(cos_psi))<=FOV
h1(ii,jj)=h1(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
end
end
end
end
end
%% Second Wall (non-RIS Wall)

for ii=1:Nx
    for jj=1:Ny
        RP=[x(ii) y(jj) -lz/2];
        % receiver position vector
        h2(ii,jj)=0;
        % reflection from North face
        for kk=1:Ny
            for ll=1:Nz
                WP2=[x(kk) -ly/2 z(ll)];
                % point of incidence in wall
                D1=sqrt(dot(TP1-WP2,TP1-WP2));
                % distance from transmitter to WP2
                cos_phi=abs(WP2(3)-TP1(3))/D1;
                cos_alpha=abs(TP1(1)-WP2(1))/D1;
                D2=sqrt(dot(WP2-RP,WP2-RP));
                % distance from WP1 to receiver
                cos_beta=abs(WP2(1)-RP(1))/D2;
                cos_psi=abs(WP2(3)-RP(3))/D2;
                if abs(acosd(cos_psi))<=FOV
                    h2(ii,jj)=h2(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end
    end
end

%% Third Wall
for ii=1:Nx
for jj=1:Ny
RP=[x(ii) y(jj) -lz/2];
% receiver position vector
h3(ii,jj)=0;
for kk=1:Ny
for ll=1:Nz
WP3=[lx/2 y(kk) z(ll)];
% point of incidence in wall
D1=sqrt(dot(TP1-WP3,TP1-WP3));
% distance from transmitter to WP1
cos_phi=abs(WP3(3)-TP1(3))/D1;
cos_alpha=abs(TP1(1)-WP3(1))/D1;
D2=sqrt(dot(WP3-RP,WP3-RP));
% distance from WP3 to receiver
cos_beta=abs(WP3(1)-RP(1))/D2;
cos_psi=abs(WP3(3)-RP(3))/D2;
if abs(acosd(cos_psi))<=FOV
h3(ii,jj)=h3(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
end
end
end
end
end

%% Fourth Wall
for ii=1:Nx
for jj=1:Ny
RP=[x(ii) y(jj) -lz/2];
% receiver position vector
h4(ii,jj)=0;
% reflection from North face
for kk=1:Ny
for ll=1:Nz
WP4=[x(kk) ly/2 z(ll)];
% point of incidence in wall
D1=sqrt(dot(TP1-WP4,TP1-WP4));
% distance from transmitter to WP4
cos_phi=abs(WP4(3)-TP1(3))/D1;
cos_alpha=abs(TP1(1)-WP4(1))/D1;
D2=sqrt(dot(WP4-RP,WP4-RP));
% distance from WP4 to receiver
cos_beta=abs(WP4(1)-RP(1))/D2;
cos_psi=abs(WP4(3)-RP(3))/D2;
if abs(acosd(cos_psi))<=FOV
h4(ii,jj)=h4(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
end
end
end
end
end
P_rec_A1=(h1+h2+h3+h4)*P_total.*Ts.*G_Con;
P_rec_A2=fliplr(P_rec_A1);
% received power from source 2, due to symmetry no need separate
% calculations
P_rec_A3=flipud(P_rec_A1);
P_rec_A4=fliplr(P_rec_A3);
P_rec_total_1ref=P_rec_A1;%+P_rec_A2+P_rec_A3+P_rec_A4;
Powers = P_rec_total_1ref;

end