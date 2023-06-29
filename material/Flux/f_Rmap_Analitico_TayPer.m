function [ct,BiotM,beta,PermK,sigmaE_new,sigmaT_new,mflu_new,velflu_new,velflu_total] = ...
    f_Rmap_Analitico_TayPer(eps_new,phi_new_dup,phi_new_n,p_new_chi,p_new_sigma,...
                                   phi_new_n1,e_DatMatSet,theta,Dtime)
   
   %Tensor constitutivo elastico
   ct = e_DatMatSet.e_DatSet.e_DatMat.ce;
   %Tensor de Biot
   BiotM = e_DatMatSet.e_DatSet.e_DatMat.m_Biot;
   %Modulo de Biot
   beta = e_DatMatSet.e_DatSet.e_DatMat.beta;
   %Tensor de permeabilidad
   PermK = e_DatMatSet.e_DatSet.e_DatMat.m_PermK;
   Omega_micro = e_DatMatSet.omegaMicro_p;
   lado_micro = Omega_micro^(0.5);
   
   In = (lado_micro^4)/12;
   Inercia_micro = [In 0 ;
                    0 In];
   %DEFICIENTE
   k=10000000;
    zita = 0;
    for j=1:k
       zita = zita + (1/((j^2)*exp(j^2)));  
    end
   %Delta de tensiones efectivas
   sigmaE_new = [ct(1,2)*eps_new(2,1); ct(2,2)*eps_new(2,1); ct(1,2)*eps_new(2,1); 0]; 
   %Tensiones totales
   sigmaT_new = sigmaE_new - BiotM*p_new_sigma;
   %Delta de contenido de masa del fluido
   mflu_new = BiotM(2,1)*eps_new(2,1) + beta*p_new_chi; 
   
   %Velocidad de filtracion del fluido por delta de tiempo al tiempo "n+theta"
   velflu_sta = [0; (-PermK(2,2)*(theta*phi_new_dup(2,1)+phi_new_n(2,1)))]; 
   omega = PermK(2,2)/beta;
   c1 = (6*zita/pi^2)*exp((-4*pi^2*omega*Dtime)/lado_micro^2);
%     c2 = beta*((lado_micro^2)/12);
   velflu_dyn = [0; ((c1-1)*(beta*Inercia_micro*phi_new_dup(2,1)/(Omega_micro*Dtime)))]; 
   velflu_new = velflu_sta + velflu_dyn;
%    c1 = (BiotM(2,1)^2)*(1/ct(2,2))*((lado_micro^2)/12);
%    c2 = beta*((lado_micro^2)/12);
%    velflu_dyn = [0; ((-c1*phi_new_dup(2,1)/Dtime)-(c2*phi_new_dup(2,1)/Dtime))]; 
%    velflu_new = velflu_sta + velflu_dyn;
   

   velflu_sta_n1 = [0; (-PermK(2,2)*phi_new_n1(2,1))]; 
   velflu_dyn_n1 = [0; ((c1-1)*(beta*Inercia_micro*phi_new_dup(2,1)/(Omega_micro*Dtime)))]; 
%    velflu_dyn_n1 = [0; ((-c1*phi_new_dup(2,1)/Dtime)-(c2*phi_new_dup(2,1)/Dtime))]; 
   velflu_total = velflu_sta_n1 + velflu_dyn_n1;

   
   

   
   