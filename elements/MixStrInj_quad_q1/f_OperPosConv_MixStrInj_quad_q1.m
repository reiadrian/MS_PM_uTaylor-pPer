function m_VarAuxElem = f_OperPosConv_MixStrInj_quad_q1(m_VarAuxElem,m_VarAuxGP,e_DatMatSet,m_CT,e_VG)

   %parfor iElem = 1:nElem
%    for iElem = 1:nElem
%       
%       condBif = m_VarAuxElem(1);
%       if ~condBif
%          %Para el an�lisis de bifurcaci�n se toma el tensor de punto de gauss 5 para el impl�cito o 10 para el
%          %caso del implex (se asume que este elemento finito siempre tiene 5 puntos de gauss).
%          if esImplex
%             %Tensores tangentes constitutivos para bifurcaci�n (se toma el tensor impl�cito).
%             m_CTBif = m_CT(:,:,10,iElem);
%          else
%             m_CTBif = m_CT(:,:,5,iElem);
%          end
%          condBif = f_CondBifct(m_CTBif,e_VG);
%       end
%       m_VarAuxElem(1) = condBif;
%       
%    end

   %Por ahora se realiza este an�lisis usando el fload.
   conshyp = e_DatMatSet.conshyp;
   switch conshyp
      case {4,5,10,12}
         siavarpg = e_DatMatSet.siavarpg;
         %Se considera que este elemento finito siempre tiene 5 puntos de gauss, y que este se utiliza para
         %determinar si el elemento completo bifurc�. El fload se considera que est� en la posici�n 2 de estos
         %modelos constitutivos.
         %Se convierte autom�ticamente a double los valores l�gicos m_VarAuxGP(4*siavarpg+2,:)>0 al ser
         %m_VarAuxElem una matriz de doubles.
         %Notar que se considera que pertenece al dominio mixto con cualquier valor positivo, y cuando hay
         %carga el�stica (0) y descarga el�stica (-1) queda excluido. Esta determinaci�n se debe realizar en
         %todos los pasos, a�n despu�s de haber detectado la bifurcaci�n ya que determinar la descarga.
         m_VarAuxElem(1,:) = m_VarAuxGP(4*siavarpg+2,:)>0;
      otherwise
         error(['Operaciones de fin de paso: Elemento finito quadril�tero Q1 mixto con inyecci�n ',...
            'de deformaci�n: Modelo constitutivo no definido'])
   end

end