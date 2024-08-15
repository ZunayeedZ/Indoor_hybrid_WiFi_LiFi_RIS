function [Powers] = Power_CAL_RIS(xx,yy,m,P_total,Adet,rho,Ts,FOV,G_Con,lx,ly,lz,XT,YT,ZT,lambda,dsep,M,N,Nx,Ny,Nz,dA,xRIS,zRIS,x,y,z)
%first transmitter calculation
TP1=[XT(1,xx,1) YT(yy,1,1) ZT(1,1,1)];
% transmitter position
%TPV=[0 0 -1];
% transmitter position vector
%RPV=[0 0 1];
% receiver position vector
%WPV1=[1 0 0];
% position vector for wall 1

%% First Wall
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
%% Second Wall (RIS Wall)
for ii=1:Nx
    for jj=1:Ny
        RP=[x(ii) y(jj) -lz/2];
        % receiver position vector
        %h2(ii,jj)=0;
        h_LR_LOS(ii,jj)=0;
        h_LR(ii,jj)=0;
        h_RU_LOS(ii,jj)=0;
        h_RU(ii,jj)=0;
        % reflection from North face
        for kk=1:N
            for ll=1:M
                % point of incidence in RIS surface
                WP2=[xRIS(ll) -ly/2 zRIS(kk)];
                % distance from transmitter to RIS surface
                D1=sqrt(dot(TP1-WP2,TP1-WP2));
                cos_phi=abs(WP2(3)-TP1(3))/D1;
                cos_alpha=abs(TP1(1)-WP2(1))/D1;
                % distance from RIS surface to receiver
                D2=sqrt(dot(WP2-RP,WP2-RP));
                cos_beta=abs(WP2(1)-RP(1))/D2;
                cos_psi=abs(WP2(3)-RP(3))/D2;
                phase_LR(ll) = exp((-1i*2*pi*(ll-1)*dsep.*cos_alpha)/lambda);
                phase_RU(ll) = exp((-1i*2*pi*(ll-1)*dsep.*cos_beta)/lambda);
                v = exp(1i.*acos(cos_alpha));
                phase_diag = diag(v);
                if abs(acosd(cos_psi))<=FOV
                    h_LR_LOS(ii,jj) =h_LR_LOS(ii,jj)+(m+1)*Adet.*cos_phi.^(m+1)./(2*pi.*D1.^2);
                    h_LR(ii,jj) = h_LR(ii,jj)+ h_LR_LOS(ii,jj).*(transpose(phase_LR(ll)));
                    h_RU_LOS(ii,jj) =h_RU_LOS(ii,jj)+(m+1)*Adet.*cos_beta.^(m+1)./(2*pi.*D2.^2);
                    h_RU(ii,jj) = h_RU(ii,jj)+ h_RU_LOS(ii,jj).*(transpose(phase_RU(ll)));
                    %h2(ii,jj)=h2(ii,jj)+(m+1)*Adet*rho*dA*cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    
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
%Finding the components of total channel gain of the RIS path
h_RU_hermitian = transpose(h_RU);
%Total Channel Gain LED - RIS - User
h2 = h_RU_hermitian*phase_diag*h_LR;
P_rec_A1=(h1+h2+h3+h4)*P_total.*Ts.*G_Con;
%P_rec_A2=fliplr(P_rec_A1);
% received power from source 2, due to symmetry no need separate
% calculations
%P_rec_A3=flipud(P_rec_A1);
%P_rec_A4=fliplr(P_rec_A3);
P_rec_total_1ref=P_rec_A1;%+P_rec_A2+P_rec_A3+P_rec_A4;
Powers = P_rec_total_1ref;

end