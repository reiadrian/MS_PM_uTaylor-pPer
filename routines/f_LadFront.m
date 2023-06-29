function  m_ConecFront = f_LadFront(m_ListEl,m_ElemVec,e_DatSet,e_VG)
   
   %Esta funci�n no funcionar correctamente el caso del grupo de elemento rodeen un agujero, ya que
   %va devolver los elemento lineales de la frontera interna, que no interviene en la integrales de
   %la homogineizaci�n de deformaci�n.
   %
   nSet = e_VG.nSet;
   nElList = length(m_ListEl);   
   %El procedimiento para la determinaci�n de los lados frontera se basa en buscar los elementos que
   %est�n en m_ListEl en la matriz de elementos vecinos (m_ElemVecList), considerando solo los
   %elementos que est�n dentro del dominio indicado y formado por los elementos por m_List, es decir
   %el grupo de elementos al que se busca la frontera. En las posiciones en que se encuentra
   %coincidencia dentro de m_ElemVecList significa que elementos de m_ListEl tienen lados comunes, y
   %por lo tanto lados que no son frontera (estos se indican con false en la matriz m_IndLadFront).
   m_ElemVecList = m_ElemVec(m_ListEl,:);
   m_IndLadFront = true(nElList,size(m_ElemVec,2));
   for iEl = 1:nElList
      %Se sigue la convenci�n de que NaN indica un lado frontera de la malla, y por lo tanto ac�
      %tambi�n (recordar que NaN==CualquierNro = false, y por lo tanto se mantiene como true en la
      %matriz m_IndLadFront).
      %Notar que los valores con cero de m_ElemVecList (son valores que corresponde a elementos que
      %menor cantidad de nodos que el m�ximo n�mero de nodos por elemento de la malla) quedan con
      %true, es decir es una lado frontera. Esto es incorrecto ya que en realidad es lado no existe,
      %pero como para buscar la conectividades de los lados de frontera se los ignora, la b�squeda
      %funciona correctamente.
      m_IndLadFront(m_ElemVecList==m_ListEl(iEl)) = false;
   end
   %Otra opci�n para calcular m_IndLadFront es (pero es m�s lento):
   %m_IndLadFront = ~ismember(m_ElemVecList,m_ListEl);
   %
   %Como se no se conoce el n�mero de lados por los valores con cero m_ElemVecList que marca lados
   %de frontera incorrectamente se almacena primero en una celda las conectividades de frontera y
   %despu�s se las junta en una matriz.
   c_ConecFront = cell(1,nSet);
   for iSet = 1:nSet
      m_Conec = e_DatSet(iSet).conec;
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      npe = e_DatSet(iSet).e_DatElem.npe;
      m_IndElemSet = e_DatSet(iSet).m_IndElemSet;
      %Ordenamiento de las conectividades para que los nodos est�n ordenados seg�n lo utilizado en
      %la funci�n f_ElemVecino.
      switch eltype
         case {2,4,8,10,20,32}
            % Tri�ngulo de 3 nodos (est�ndar y SDA), Cuadr�ngulo de 4 nodos (est�ndar y bbar).
            m_NodLadElem = m_Conec;
         otherwise
            error(['Determinaci�n de lados frontera de un grupo de elementos: Elemento finito ',...
               'no definido.'])
      end
      m_IndListElSet = bsxfun(@eq,m_ListEl',m_IndElemSet);
      [m_IndLadFrontSetFila,m_IndLadFrontSetCol] = find(m_IndLadFront(any(m_IndListElSet,2),1:npe));
      m_NodLadElem = m_NodLadElem(any(m_IndListElSet,1),:);
      nLadFront = length(m_IndLadFrontSetFila);
      m_ConecFront = zeros(2,nLadFront);
      for iLadFront = 1:nLadFront
         iFil = m_IndLadFrontSetFila(iLadFront);
         iCol = m_IndLadFrontSetCol(iLadFront);
         m_ConecFront(1,iLadFront) = m_NodLadElem(iFil,iCol);
         iCol = iCol+1;
         iCol(iCol>npe) = 1;
         m_ConecFront(2,iLadFront) = m_NodLadElem(iFil,iCol);
      end
      c_ConecFront{iSet} = m_ConecFront;
   end
   m_ConecFront = [c_ConecFront{:}]';

end