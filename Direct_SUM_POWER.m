function [direct_power] = Direct_SUM_POWER(xxd,yyd,m,P_total,Adet,Ts,FOV,G_Con,lz,XTd,YTd,XRd,YRd)

D1d=sqrt((XRd-XTd(1,xxd)).^2+(YRd-YTd(yyd,1)).^2+lz^2);
% distance vector from source 1
cosphi_A1d=lz./D1d;
% angle vector
receiver_angle=acosd(cosphi_A1d);
H_A1d=(m+1)*Adet.*cosphi_A1d.^(m+1)./(2*pi.*D1d.^2);
% channel DC gain for source 1
P_rec_A1d=P_total.*H_A1d.*Ts.*G_Con;
% received power from source 1;
P_rec_A1d(find(abs(receiver_angle)>FOV))=0;
% if the anlge of arrival is greater than FOV, no current is generated at the photodiode.
%P_rec_A2d=fliplr(P_rec_A1d);
% received power from source 2, due to symmetry no need separate calculations
%P_rec_A3d=flipud(P_rec_A1d);
%P_rec_A4d=fliplr(P_rec_A3d);
direct_power=P_rec_A1d;

end