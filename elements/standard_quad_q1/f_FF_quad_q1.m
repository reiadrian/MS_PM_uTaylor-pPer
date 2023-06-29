function m_FFe = f_FF_quad_q1(e_DatElemSet,e_VG)
   
   ndn = e_VG.ndn;
   struhyp = e_VG.struhyp;
   dofpe = e_DatElemSet.dofpe;
   nPG = e_DatElemSet.npg;
   xg = e_DatElemSet.xg;
  
   E = xg(:,1);
   n = xg(:,2);

% Funciones de forma
   m_FFe = zeros(ndn,dofpe,nPG);
%AA: Agrege if y FF para Elemento de 8 nodos y además 24 gdl (tomar en cuenta la poropresion)
   
if dofpe==8  %Elemento de 4 nodos
        
% Funciones de forma: Elemento cuadrilatero de 4 nodos
        N1 = 1/4*(1-E).*(1-n);
        N2 = 1/4*(1+E).*(1-n);
        N3 = 1/4*(1+E).*(1+n);
        N4 = 1/4*(1-E).*(1+n);
   
   switch struhyp
      case {1,2,20} %Deformación y tensión plana, y Large Deformations with plane deformation.
         %m_FFpg = [N1,0,N2,0,N3,0,N4,0;0,N1,0,N2,0,N3,0,N4];
         m_FFe(1,1,:) = N1;
         m_FFe(1,3,:) = N2;
         m_FFe(1,5,:) = N3;
         m_FFe(1,7,:) = N4;
         m_FFe(2,2,:) = N1;
         m_FFe(2,4,:) = N2;
         m_FFe(2,6,:) = N3;
         m_FFe(2,8,:) = N4; %AA: Cambie N3 por N4
       otherwise %struhyp
         error('Elemento Cuad_Q1: Matriz de forma: Hipótesis estructural no definido.');
   end %struhyp
   
elseif dofpe==16  %AA: Agrege FF para Elemento de 8 nodos. Problema monofase solida
    
% Funciones de forma: Elemento cuadrilatero de 8 nodos
        N1 = 1/4*(1-E).*(1-n).*(-E-n-1);
        N2 = 1/4*(1+E).*(1-n).*(E-n-1) ;
        N3 = 1/4*(1+E).*(1+n).*(E+n-1) ;
        N4 = 1/4*(1-E).*(1+n).*(-E+n-1);
        N5 = 1/2*(1-E.^2).*(1-n)       ;
        N6 = 1/2*(1+E).*(1-n.^2)       ;
        N7 = 1/2*(1-E.^2).*(1+n)       ;
        N8 = 1/4*(1-E).*(1-n.^2)       ;

   switch struhyp
      case {1,2,20} %Deformación y tensión plana, y Large Deformations with plane deformation.
        m_FFe(1,1,:)  = N1;
        m_FFe(1,3,:)  = N2;
        m_FFe(1,5,:)  = N3;
        m_FFe(1,7,:)  = N4;
        m_FFe(1,9,:)  = N5;
        m_FFe(1,11,:) = N6;
        m_FFe(1,13,:) = N7;
        m_FFe(1,15,:) = N8;
        
        m_FFe(2,2,:)  = N1;
        m_FFe(2,4,:)  = N2;
        m_FFe(2,6,:)  = N3;
        m_FFe(2,8,:)  = N4; 
        m_FFe(2,10,:) = N5;
        m_FFe(2,12,:) = N6;
        m_FFe(2,14,:) = N7;
        m_FFe(2,16,:) = N8;
       otherwise %struhyp
        error('Elemento Cuad_Q1: Matriz de forma: Hipótesis estructural no definido.');
   end %struhyp
   
% elseif dofpe==24  %AA: Agrege FF para Elemento de 8 nodos. Problema bifase (Suelo saturado)
%                   %    Incluye 3 gdl (2 desplazamientos y 1 poropresion)
%     
% % Funciones de forma: Elemento cuadrilatero de 4 nodos (P/poropresiones)
%         SF1 = 1/4*(1-E).*(1-n);
%         SF2 = 1/4*(1+E).*(1-n);
%         SF3 = 1/4*(1+E).*(1+n);
%         SF4 = 1/4*(1-E).*(1+n);
%         
% % Funciones de forma: Elemento cuadrilatero de 8 nodos (P/desplazamientos)
%         N1 = 1/4*(1-E).*(1-n).*(-E-n-1);
%         N2 = 1/4*(1+E).*(1-n).*(E-n-1) ;
%         N3 = 1/4*(1+E).*(1+n).*(E+n-1) ;
%         N4 = 1/4*(1-E).*(1+n).*(-E+n-1);
%         N5 = 1/2*(1-E.^2).*(1-n)       ;
%         N6 = 1/2*(1+E).*(1-n.^2)       ;
%         N7 = 1/2*(1-E.^2).*(1+n)       ;
%         N8 = 1/4*(1-E).*(1-n.^2)       ;
%         
%    switch struhyp
%       case {1,2,20} %Deformación y tensión plana, y Large Deformations with plane deformation.
% 
% % Desplazamientos       
%         m_FFe(1,1,:)   = N1; 
%         m_FFe(1,4,:)   = N2;
%         m_FFe(1,7,:)   = N3;
%         m_FFe(1,10,:)  = N4;
%         m_FFe(1,13,:)  = N5;
%         m_FFe(1,16,:)  = N6;
%         m_FFe(1,19,:)  = N7;
%         m_FFe(1,22,:)  = N8;
%         
%         m_FFe(2,2,:)   = N1;
%         m_FFe(2,5,:)   = N2;
%         m_FFe(2,8,:)   = N3;
%         m_FFe(2,11,:)  = N4; 
%         m_FFe(2,14,:)  = N5;
%         m_FFe(2,17,:)  = N6;
%         m_FFe(2,20,:)  = N7;
%         m_FFe(2,23,:)  = N8;
%         
% %Poropresiones
%         m_FFe(3,3,:)   = SF1; 
%         m_FFe(3,6,:)   = SF2;
%         m_FFe(3,9,:)   = SF3;
%         m_FFe(3,12,:)  = SF4; 
%         m_FFe(3,15,:)  = 0;
%         m_FFe(3,18,:)  = 0;
%         m_FFe(3,21,:)  = 0;
%         m_FFe(3,24,:)  = 0;
% 
%       otherwise
%         error('Elemento Cuad_Q1: Matriz de forma: Hipótesis estructural no definido.');
%    end
else %dofpe
    error('Elemento Cuad_Q1: Matriz de forma: Tipo de elemento no definido.');
    
end %if: dofpe %AA: Agrege FF para Elemento de 8 nodos
