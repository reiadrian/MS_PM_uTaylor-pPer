function [m_FF_d,m_FF_p] = f_FF_Bifase_quadbq(e_DatElemSet,e_VG)
   
   %ndn = e_VG.ndn;
   struhyp = e_VG.struhyp;
   %npe = e_DatElemSet.npe;
   nPG = e_DatElemSet.npg;
   xg = e_DatElemSet.xg;
   
   ndn_d = e_DatElemSet.ndn_d;
   ndn_p = e_DatElemSet.ndn_p;
   dofpe_d = e_DatElemSet.dofpe_d;
   dofpe_p=  e_DatElemSet.dofpe_p;
   
   E = xg(:,1);
   n = xg(:,2);

% Funciones de forma
%    ndn_d = ndn-1;
%    ndn_p = 1;
%    dofpe_d = ndn_d*npe;
%    dofpe_p = ndn_p*(npe-4); %AA: ver si se puede hacer general
   m_FF_d = zeros(ndn_d,dofpe_d,nPG); % Desplazamientos
   m_FF_p = zeros(ndn_p,dofpe_p,nPG); % Poropresiones 
% Funciones de forma: Elemento cuadrilatero de 4 nodos (P/poropresiones)
   SF1 = 1/4*(1-E).*(1-n);
   SF2 = 1/4*(1+E).*(1-n);
   SF3 = 1/4*(1+E).*(1+n);
   SF4 = 1/4*(1-E).*(1+n);
        
% Funciones de forma: Elemento cuadrilatero de 8 nodos (P/desplazamientos)
   N1 = 1/4*(1-E).*(1-n).*(-E-n-1);
   N2 = 1/4*(1+E).*(1-n).*(E-n-1) ;
   N3 = 1/4*(1+E).*(1+n).*(E+n-1) ;
   N4 = 1/4*(1-E).*(1+n).*(-E+n-1);
   N5 = 1/2*(1-E.^2).*(1-n)       ;
   N6 = 1/2*(1+E).*(1-n.^2)       ;
   N7 = 1/2*(1-E.^2).*(1+n)       ;
   N8 = 1/4*(1-E).*(1-n.^2)       ;
        
   switch struhyp
        case {1,2,20} %Deformacion y tension plana, y Large Deformations with plane deformation.

% P/Desplazamientos       
            m_FF_d(1,1,:)  = N1;
            m_FF_d(1,3,:)  = N2;
            m_FF_d(1,5,:)  = N3;
            m_FF_d(1,7,:)  = N4;
            m_FF_d(1,9,:)  = N5;
            m_FF_d(1,11,:) = N6;
            m_FF_d(1,13,:) = N7;
            m_FF_d(1,15,:) = N8;

            m_FF_d(2,2,:)  = N1;
            m_FF_d(2,4,:)  = N2;
            m_FF_d(2,6,:)  = N3;
            m_FF_d(2,8,:)  = N4; 
            m_FF_d(2,10,:) = N5;
            m_FF_d(2,12,:) = N6;
            m_FF_d(2,14,:) = N7;
            m_FF_d(2,16,:) = N8;

    % P/Poropresiones
            m_FF_p(1,1,:)  = SF1; 
            m_FF_p(1,2,:)  = SF2;
            m_FF_p(1,3,:)  = SF3;
            m_FF_p(1,4,:)  = SF4; 
        otherwise %struhyp
            error('Elemento Cuad_Q1: Matriz de forma: Hipotesis estructural no definido.');
   end %struhyp
%if: dofpe %AA: Agrege FF para Elemento de 8 nodos
