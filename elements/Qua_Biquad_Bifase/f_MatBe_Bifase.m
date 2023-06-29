function [m_Be_d,m_DetJe_d,m_DerCae_p,m_DetJe_p] = f_MatBe_Bifase(coord_n,e_DatElemSet,e_VG)

   %xg = e_VG.xg;
   ntens = e_VG.ntens;
   ndime = e_VG.ndime;
   %ndn = e_VG.ndn;
   struhyp = e_VG.struhyp;
   
   dofpe_d = e_DatElemSet.dofpe_d;
   dofpe_p=  e_DatElemSet.dofpe_p;
   xg = e_DatElemSet.xg;   
   npg = e_DatElemSet.npg;

% Desplazamientos
    m_Be_d = zeros(ntens,dofpe_d,npg);
    m_DetJe_d = zeros(npg,1);
    
% Poropresi�n
    m_DerCae_p = zeros(ndime,dofpe_p,npg);
    m_DetJe_p = zeros(npg,1);
   
   for iPG = 1:npg
      switch struhyp
         case {1,2}
            [m_Be_d(:,:,iPG),m_DetJe_d(iPG),m_DerCae_p(:,:,iPG),m_DetJe_p(iPG)]...
               = matrixB_Bifase_epd(coord_n,xg(iPG,:),e_DatElemSet,e_VG); 
         case 3
            fprintf('Estado 3D no disponible');
         case 4
            fprintf('Estado axisimetrico no disponible');
            %[m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_quad_q1_axm(coord_n,xg(iPG,:),e_DatElemSet,e_VG);
         case 20   %Large deformations with plane deformation.
            fprintf('Estado EPD Grandes deformaciones no disponible');
            %[m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_quad_q1_epd_ld(coord_n,xg(iPG,:),e_DatElemSet,e_VG);
         otherwise
            error('Elemento quad_q1: Matriz de deformacion B Global: Estado no disponible.');
      end
   end
   
end