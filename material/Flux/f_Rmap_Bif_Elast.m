function [ct,BiotM,beta,PermK,sigmaE_new,sigmaT_new] = ...
    f_Rmap_Bif_Elast(eps_new,e_DatMatSet,up,N4)
   
   ct = e_DatMatSet.ce;
   BiotM = e_DatMatSet.m_Biot;
   beta = e_DatMatSet.beta;
   PermK = e_DatMatSet.m_PermK;
   
   sigmaE_new = ct*eps_new; %Tensiones efectivas
   sigmaT_new = sigmaE_new - BiotM*N4*up; % Tensiones totales
   
   

   
   