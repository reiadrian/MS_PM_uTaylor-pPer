function m_ElemVec = f_ElemVecino(e_DatSet,e_VG)

   nElem = e_VG.nElem;
   nSet = e_VG.nSet;   
   
   %Se hace un bucle sobre los sets para rearmar una matriz conectividad completa. Estas matrices
   %son necesarios para buscar los lados vecinos con diferentes elementos, ya que permite una
   %b�squeda en toda la malla de los lados vecinos.
   %Ver si no traer como una variable m_npeMax con e_VG.
   npeMax = max(arrayfun(@(x)x.e_DatElem.npe,e_DatSet));
   m_NodLadElem = zeros(nElem,npeMax);
   m_NumIntNod = zeros(nElem,npeMax+2);
   for iSet = 1:nSet
      conec = e_DatSet(iSet).conec;
      nElemSet = e_DatSet(iSet).nElem;
      m_IndElemSet = e_DatSet(iSet).m_IndElemSet;
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      npe = e_DatSet(iSet).e_DatElem.npe;
      %Para considerar que la numeraci�n interna de los nodos de los elementos no sigue un orden de
      %rotaci�n (horario o antihorario) se utiliza un paso previo para reordenar seg�n los los lados
      %del elemento. 
      %Tambi�n podr�a servir si se usa la malla compatibles (continuidad de desplazamiento entre
      %elementos), para buscar por menos nodos si es una lado es vecino a otro.
      %Si la malla no fuera compatible habr�a que analizar como realizarlo (mezclar por ejemplo
      %elementos 2D cuadr�ticos y lineales).
      %Tambi�n se usa para crear una matriz que tiene la dimensi�n del n�mero de nodos por elemento
      %(npe) m�ximo de toda la malla. Si npe es menor al npeMax, en las posiciones finales se deja
      %cero.
      switch eltype
         case {2,4,8,10,20,32}
            % Tri�ngulo de 3 nodos (est�ndar y SDA), Cuadr�ngulo de 4 nodos (est�ndar, bbar y mixto con iny.
            %de deformac�n).
            m_NodLadElem(m_IndElemSet,1:npe) = conec;
         case 7
            error(['Determinac�n de elementos vecinos: No implementado para elementos 3D, solo ',...
               'para elementos 2D'])
         otherwise
            error('Determinaci�n de elementos vecinos: Elemento finito no definido.')
      end
      % Matriz de numeraci�n interna de los nodos del elemento)
      %Se repite el n�mero interno del nodo �ltimo al principio y el �ltimo al final para facilitar
      %la indexaci�n.
      m_NumIntNod(m_IndElemSet,1:npe+2) = repmat([npe,1:npe,1],nElemSet,1);
      %Me parece que todas las operaciones se podr�a realizar creando dos matrices, una indicando el
      %elemento que corresponde con la fila de otra matriz de 3 columnas. En la primera columna de
      %esta matriz se indica el n�mero de nodo previo, la segunda el nodo actual, y en la tercera el
      %siguiente.
      %Entonces el procedimiento es buscar para un cierto elemento y nodo donde tal que se verifique
      %el nodo actual y el siguiente sea igual al nodo inicial y al actual (respectivamente).
      %Esto evitar�a crear matrices con ceros para elementos cuadr�ngulares y triangulares lineales.
      %Ver si as� ser�a m�s r�pido.
   end
   %
   %Se define que cuando un elemento no tiene en alg�n lado un vecino, lleva un NaN en ese lado.
   %Cuando el elemento tiene menos nodos por elementos que npeMax, coloca un cero
   m_ElemVec = nan(nElem,npeMax);
   %Como cuando se determina el elemento vecino de un elemento, tambi�n se sabe un vecino del 
   %primero, en cada lado que se analiza se guarda en las filas de los dos elementos el
   %correspondiente vecino. Para no pasar dos veces por el mismo, se guarda una lista de los
   %elementos y lados analizados
   m_Analiz = false(nElem,npeMax);
   %Si se usara parfor habr�a que sacar unas optimizaciones de no buscar el mismo lados 2 veces.
   for iElem = 1:nElem
      %Notar que el primer valor de m_NumIntNod es igual al n�mero de nodos del elemento.
      npe = m_NumIntNod(iElem,1);
      for iNodInLad = 1:npe
         % Se salta los lados analizados.
         if ~m_Analiz(iElem,iNodInLad)
            % Encuentra los elementos y el lado que tienen el nodo inicial del lado analizado.
            nodInLad = m_NodLadElem(iElem,iNodInLad);
            [m_ElemVecNodFiLad,m_iNodFiLadVec] = find(m_NodLadElem==nodInLad);
            % �ndice del nodo final del lado que se est� analizando.
            %+1 por que se quiere tomar el nodo siguiente.
            %iNodEl+1+1 (+1, por la diferencia entre m_NodLadElem y m_NumIntNod)
            iNodFiLad = m_NumIntNod(iElem,iNodInLad+2);
            nodFiLad = m_NodLadElem(iElem,iNodFiLad);
            %Este loop capaz se puede vectorizar, ver si vale la pena.
            for iNodFi = 1:length(m_ElemVecNodFiLad)
               %Ac� se podr�a descartar el caso de elVec==iElem, pero se evitar�a pocas operaciones
               %y agregar�a un condicional por lo que posiblemente no vale la pena.
               elVec = m_ElemVecNodFiLad(iNodFi);
               % N�mero de �ndice de los nodos global final de los lados que tienen igual nodo inicial.
               %Se utiliza -1 ya que el lado com�n de dos elementos tienen sentidos opuestos de
               %numeraci�n local, por lo que se debe buscar el nodo inicial del lado.
               %m_IndNodInLad(iNodIn)+1-1 (+1, por la diferencia entre m_NodLadElem y m_NumIntNod)
               iNodInLadVec = m_NumIntNod(elVec,m_iNodFiLadVec(iNodFi));  
               nodInLadVec = m_NodLadElem(elVec,iNodInLadVec);
               % Determinaci�n de elementos que tienen el lado en com�n.
               %Se ve si tiene un lado contiguo.
               if nodFiLad==nodInLadVec
                  %Se almacena el vecino.
                  m_ElemVec(iElem,iNodInLad) = elVec;
                  %Se almacena el elemento analizado en el vecino.
                  m_ElemVec(elVec,iNodInLadVec) = iElem;
                  %Se guarda lado y elemento analizado, para no volver analizar el vecino (un
                  %elemento solo puede tener un elemento vecino por cada lado).
                  m_Analiz(iElem,iNodInLad) = true;
                  m_Analiz(elVec,iNodInLadVec) = true;
                  break
               end
            end
         end
         %Esto no ser�a necesario ya que no vuelve analizarse el lado de este elemento
         %nuevamente (es decir no se pasa de nuevo por la posici�n (iElem,iNodEl)).
         m_Analiz(iElem,iNodInLad) = true;
      end
   end
   
end
            


      