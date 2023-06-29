function e_DatSet = f_MatBT(xx,e_DatSet,e_VG)
   
   ntens = e_VG.ntens;
   ndime = e_VG.ndime;
   nSet = e_VG.nSet;
   protype= e_VG.protype;
   %struhyp = e_VG.struhyp;
   
   for iSet = 1:nSet
      
      conec = e_DatSet(iSet).conec;
      nElem = e_DatSet(iSet).nElem;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      eltype = e_DatElemSet.eltype;
      npe = e_DatElemSet.npe;
      dofpe = e_DatElemSet.dofpe;
      npg = e_DatElemSet.npg;
      wg = e_DatElemSet.wg;
      
      %La evaluacion de las funciones de forma en los puntos de gauss no depende de la forma del elemento, por
      %lo que se evalua para todos los puntos de gauss del elemento master de todos los elementos, ahorrando
      %memoria.
      switch protype
          case 0
             switch eltype
                 case {2,10,32}
                    m_FF = f_FF_tria_t1(e_DatElemSet,e_VG);
                 case {4,8,16,20,21,22,23,31,108} %AA: 16
                    m_FF = f_FF_quad_q1(e_DatElemSet,e_VG);
                 otherwise
                   error('Matrices de funcion de forma: Elemento no implementado.')
             end %eltype
          case 1
% AA: creo funcion para considerar en ella la posibilidad de dos FF
% Una de 8 nodos p/desp y 4 nodos para poropresion (Independiza de eltype,
% ya que deberia entrar para 16 solamente (Ya que en caso de suelo saturado, 
% minimamente ya considero el elemento de 8 nodos). BUSCAR SOLUCION MAS CONVENIENTE
             switch eltype
                 case 16
                    [m_FF_d,m_FF_p] = f_FF_Bifase_quadbq(e_DatElemSet,e_VG); %AA: cree funcion
                 otherwise
                    error('Matrices de funcion de forma: Elemento no implementado.') 
             end  %eltype
          otherwise %protype
             error('Lectura de datos: Hipotesis de tipo de problema no implementada')
      end %protype
      %
      switch protype %AA
          case 0
              m_BT = zeros(ntens,dofpe,npg,nElem);
              m_DetJT = zeros(npg,nElem);
          case 1 
              dofpe_d = e_DatElemSet.dofpe_d;
              dofpe_p=  e_DatElemSet.dofpe_p;
              m_BT_d = zeros(ntens,dofpe_d,npg,nElem);
              m_DetJT_d = zeros(npg,nElem);
              m_DerCa_p= zeros(ndime,dofpe_p,npg,nElem); % Cuidado con ndime
              m_DetJT_p = zeros(npg,nElem);
          otherwise %protype
              error('Lectura de datos: Hipotesis de tipo de problema no implementada')   
      end %protype
      
      if eltype==31 
         dN_xy = zeros(ndime,npe,nElem);
      end
      
      %parfor iElem = 1:nElem
      for iElem = 1:nElem
         %Para que funcione correctamente con los indices, el parfor se hace en una funcion separada 
         %que resuelva los loops sobre los puntos de Gauss

         coord_n = f_CoordElem(xx,conec(iElem,:));         
         switch protype %AA
             case 0 %AA
                 switch eltype
                    case {2,10}
                       [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_tria_t1(...
                          coord_n,e_DatElemSet,e_VG);
                    case {4,16,20,21,22,23,108} %AA: 16
                      %####################################################################################
                       [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_quad_q1(...
                          coord_n,e_DatElemSet,e_VG);
                      %####################################################################################
                    case 5
                       [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_barra2D(coord_n,e_DatElemSet,e_VG);
                    case 7
                       error('Matrices de deformacion: Falta implementar el caso de hexaedros de 8 nodos')
                    case 8
                       [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_bbar_q1(coord_n,e_DatElemSet,e_VG);
                    case 31
                       [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = ...
                           f_MatBe_quad_q1(coord_n,e_DatElemSet,e_VG);
                       E = 0; n =0;
                       dN_E = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
                       dN_n = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
                       dN_En = [dN_E ; dN_n];
                       J11 = coord_n(1,:)*dN_E';
                       J12 = coord_n(2,:)*dN_E';
                       J21 = coord_n(1,:)*dN_n';
                       J22 = coord_n(2,:)*dN_n';
                       J = [J11 J12 ; J21 J22];
                       dN_xy(:,:,iElem) = J\dN_En;
                    case 32
                       [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_tria_t1(coord_n,e_DatElemSet,e_VG);
                       dN_E = [-1 1 0];
                       dN_n = [-1 0 1];
                       dN_En = [dN_E ; dN_n];
                       J11 = coord_n(1,:)*dN_E';
                       J12 = coord_n(2,:)*dN_E';
                       J21 = coord_n(1,:)*dN_n';
                       J22 = coord_n(2,:)*dN_n';
                       J = [J11 J12 ; J21 J22];
                       dN_xy(:,:,iElem) = J\dN_En;
                     otherwise %eltype
                       error('Matrices de deformacion: Elemento no implementado.')
                 end %eltype
             case 1 %AA
                switch eltype
                    case 16
                        %####################################################################################
                       [m_BT_d(:,:,:,iElem),m_DetJT_d(:,iElem),m_DerCa_p(:,:,:,iElem),m_DetJT_p(:,iElem)]...
                           = f_MatBe_Bifase(coord_n,e_DatElemSet,e_VG); %AA: cree funcion
                       %####################################################################################
                    otherwise %eltype
                       error('Matrices de funcion de forma: Elemento no implementado.') 
                end  %eltype
             otherwise %protype
                error('Lectura de datos: Hipotesis de tipo de problema no implementada')
          end %protype
      end %for
      %
      switch protype %AA
          case 0
             e_DatSet(iSet).m_BT = m_BT;
             e_DatSet(iSet).m_DetJT = m_DetJT;
             e_DatSet(iSet).m_FF = m_FF;
          case 1
             e_DatSet(iSet).m_BT_d = m_BT_d;
             e_DatSet(iSet).m_DetJT_d = m_DetJT_d;
             e_DatSet(iSet).m_DerCa_p = m_DerCa_p;
             e_DatSet(iSet).m_DetJT_p = m_DetJT_p;
             e_DatSet(iSet).m_FF_d = m_FF_d;
             e_DatSet(iSet).m_FF_p = m_FF_p;
          otherwise %protypes
             error('Matrices no implementadas: Elemento no implementado.')
      end % protypes
      %Volumen de los elementos
      %El punto de gauss ya viene multiplicado por el espesor.
      %Para considerar que el elemento SDA_tria_t1 tiene dos PG ubicados en la misma posicion para
      %calculo del volumen del elemento se debe considerar solo de ellos (ver si no poner que el
      %punto de gauss singular tenga peso nulo).
      %Parecido ocurre con el elemento MixStrInj_quad_q1, donde el PG 5 tiene el peso como si fuera un unico
      %PG del elemento, por lo que para determinar el volumen se utiliza los 4 primeros PGs.
      switch eltype %(Con 16 dirige a caso Bifase 0 y 1)
          case {2,4,7,8,16,108} %AA
              switch protype %AA
                   case 0
                       e_DatSet(iSet).m_VolElem = wg'*m_DetJT;
                   case 1
                       e_DatSet(iSet).m_VolElem_d = wg'*m_DetJT_d; % VER CUESTION DE ESPESOR YA IMPUESTA
                       e_DatSet(iSet).m_VolElem_p = wg'*m_DetJT_p;
              end %protype
          case 10
              e_DatSet(iSet).m_VolElem = wg(1)'*m_DetJT(1,:);
           case {20,21,22,23}
              e_DatSet(iSet).m_VolElem = wg(1:4)'*m_DetJT(1:4,:);
          case {31,32}
              e_DatSet(iSet).m_VolElem = wg'*m_DetJT;
              e_DatSet(iSet).dN_xy = dN_xy;
          otherwise %eltype
            error('Volumen del elemento: Elemento no implementado.')
      end %eltype
      
   end %for(iSet)

end
