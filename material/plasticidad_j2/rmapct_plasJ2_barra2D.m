function [ct] = rmapct_plasJ2_barra2D (eps_new,hvar_old,aux_var,Eprop,ce,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING MATRIZ TANGENTE PARA BARRAS EN 2D                                     *
%*  MODELO DE PLASTICIDAD J2 CON ENDURECIMIENTO ISOTROPO LINEAL                           *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global ntens
ntens = e_VG.ntens;

% Variables
% ���������
eps_old_plast  = hvar_old(1:ntens,1);
alpha_old      = hvar_old(ntens+1,1);
E              = Eprop(4);
sigmay         = Eprop(6);
tita           = Eprop(8);
hb             = Eprop(9);
kinfb          = Eprop(10);
kcerob         = Eprop(11);
delta          = Eprop(12);

% C�lculo del estado trial
% ������������������������
eps_trial_plast = eps_old_plast;
alpha_trial     = alpha_old;
sigma_trial     = ce*(eps_new - eps_old_plast);

% Evaluaci�n del criterio de fluencia
% �����������������������������������
k         = hb;
f_trial   = abs(sigma_trial) - (sigmay + k*alpha_trial);

if f_trial <= 0
   
   % Paso el�stico
   % �������������
   ct = ce;
   
else
   
   % Paso pl�stico
   % �������������
   ct = E*k/(E+k);
   
end