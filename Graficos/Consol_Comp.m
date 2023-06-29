clear;
clc;
close all
% Tv = 0.1:0.1:1;
% T = 0.10496*10^-2;
T = 0.10496;
% T = 0.10496*10^-4;
Tv = [0.088 0.34 0.76];
z_H = 0:0.05:1;
m = 0:100000;
k1 = size(Tv,2);
k2 = size(z_H,2);
k3 = size(m,2);

% for ik1=1:k1
%     p_p0 = zeros(k2,1);
%     for ik2=1:k2
%         for ik3=1:k3
%             M = 0.5*pi()*(2*m(ik3)+1);
%             Mz_H = M*z_H(ik2);
%             sMz_H = sin(Mz_H);
%             eTv = exp(1)^(-(M^2)*Tv(ik1));
%             p_p0(ik2,1) = p_p0(ik2,1) + ((2/M)*sMz_H*eTv);
%         end
%     end
%     figure(1)
%     plot(p_p0,z_H/2)
%     hold on
% end
% 
% grid on
% grid minor
% ax = gca;
% ax.YDir = 'reverse';
A = Tv/T;
cv1 = Tv(1)*(7^2)/A(1);
cv2 = Tv(2)*(7^2)/A(2);
% cv3 = Tv(3)*(7^2)/A(3);
gamma_w=9.81;
E=287281.6/1000;
nu=0.4;
e0=0.9;
k=1.2192e-6;
% gamma_w=62.43;
% E=6000;
% nu=0.3;
% e0=0.9;
% k=4e-6;
factor=1-((2*(nu^2))/(1-nu));
Em = E/factor;
K=(3*E)/(3*(1-2*nu));
k_s = k/86400;
mv=(1/Em);%/(1+e0);
cvd=k/(mv*gamma_w);
cv_H = cvd/((7*0.3048)^(2));
% cv1_2 = Tv(1)*(7*0.3048)^(2);
% cv2_2 = Tv(2)*(7*0.3048)^(2);
% cv3_2 = Tv(3)*(7*0.3048)^(2);