function [ct,BiotM,beta,PermK,sigmaE_new,sigmaT_new,mflu_new,velflu_sta,mflu_n1,velflu_n1]= ...
    f_Rmap_Bif_Elast_Monoscale(eps_new,phi_new_dup,phi_new_n,p_new_chi,p_new_sigma,...
                                   eps_n1,phi_new_n1,p_n1,e_DatMatSet,theta)
   

   %Tensor constitutivo elastico
   ct = e_DatMatSet.ce;
   %Tensor de Biot
   BiotM = e_DatMatSet.m_Biot;
   %Modulo de Biot
   beta = e_DatMatSet.beta;
   %Tensor de permeabilidad
   PermK = e_DatMatSet.m_PermK;
      
   %Delta de tensiones efectivas
   sigmaE_new = ct*eps_new; 
   %Delta de tensiones totales
   sigmaT_new = sigmaE_new - BiotM*p_new_sigma; 
   %Delta de contenido de masa del fluido
   mflu_new = BiotM'*eps_new + beta*p_new_chi;
   %Velocidad de filtracion del fluido por delta de tiempo al tiempo "n+theta"
   velflu_sta = -PermK*(theta*phi_new_dup+phi_new_n); 
   
      %Contenido de masa del fluido
   mflu_n1 = BiotM'*eps_n1 + beta*p_n1; 
   %Velocidad de filtracion del fluido por delta de tiempo al tiempo "n+theta"
   velflu_n1 = -PermK*phi_new_n1; 

%    velflu_new = -PermK*phi_new; %VER PORQUE NO ES MENOS (-)
%    velflu_new = -PermK*DerivN*phi_new; % Velocidad de filtracion del fluido
%    mflu_new = BiotM'*M_Id*B*ud;

   


   
   

   
   