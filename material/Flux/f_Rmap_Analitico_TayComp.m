function [ct,BiotM,beta,PermK,sigmaE_new,sigmaT_new,mflu_new,velflu_new,velflu_total] = ...
    f_Rmap_Analitico_TayComp(eps_new,phi_new_dup,phi_new_n,p_new_chi,p_new_sigma,...
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
   
   %Delta de tensiones efectivas
   sigmaE_new = ct*eps_new; 
   %Tensiones totales
   sigmaT_new = sigmaE_new - BiotM*p_new_sigma;
   %Delta de contenido de masa del fluido
   mflu_new = BiotM'*eps_new + beta*p_new_chi; 
   
   %Velocidad de filtracion del fluido por delta de tiempo al tiempo "n+theta"
   velflu_sta = -PermK*(theta*phi_new_dup+phi_new_n); 
   velflu_dyn = (beta*Inercia_micro*phi_new_dup)/(Dtime*Omega_micro); 
   velflu_new = velflu_sta - velflu_dyn;
   

   velflu_sta_n1 = -PermK*phi_new_n1; 
   velflu_dyn_n1 = (beta*Inercia_micro*phi_new_dup)/(Dtime*Omega_micro); 
   velflu_total = velflu_sta_n1 - velflu_dyn_n1;


   
   

   
   