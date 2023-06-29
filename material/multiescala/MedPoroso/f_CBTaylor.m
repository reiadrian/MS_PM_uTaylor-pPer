function m_CondBord = f_CBTaylor(m_ElemLoc,e_DatSet,m_CondBord,doff,e_VG)

   ndn = e_VG.ndn;
   nNod = e_VG.nnod;
   %
   %Grado de libertad con dezplazamiento prescripto.
   nroGdl = 1;
   %
   m_NodDomLoc = unique(cell2mat(arrayfun(@(x)reshape(x.conec(m_ElemLoc(x.m_IndElemSet),:),[],1),e_DatSet,...
      'UniformOutput',false)))';
   nNodDomLoc = length(m_NodDomLoc);
   m_DoffNod = reshape(doff,ndn,nNod);
   if any(any(m_DoffNod(m_NodDomLoc)))
      %Las condiciones de borde van a pisar las previas, excepto si es una condici�n de m�nima restricci�n
      %externa (tipo 3). Por eso se lee qu� nodos y gdl tienen condiciones impuestas.
      %Se transforma la indicaci�n �nica de grados de libertad restringidos en columnas separadas.
      %Se ordena en la filas ndn y en las columnas los nodos con restricci�n.
      m_DirRestr = zeros(ndn,size(m_CondBord,1));
      m_DirRestr(ndn,:) = m_CondBord(:,2);
      for iGdl = 1:ndn-1
         decGdl = 10^(ndn-iGdl);
         m_DirRestrGdl = fix(m_DirRestr(ndn,:)/decGdl);
         m_DirRestr(iGdl,:) = m_DirRestrGdl;
         m_DirRestr(ndn,:) = m_DirRestr(ndn,:)-m_DirRestrGdl*decGdl;
      end
      m_NodRestr = m_CondBord(:,1)';
      %Se detecta los grados de libertad que tiene m�nima restricci�n tipo 3
      m_DoffNodMinRest = false(ndn,nNod);
      m_DoffNodMinRest(:,m_NodRestr) = m_DirRestr==3;
      %Se fija los grados de libertad que deber�a condici�n con desplazamiento prescripto por la condici�n de
      %borde de Taylor en el dominio micro localizado.
      m_DoffNodTaylor = false(ndn,nNod);
      m_DoffNodTaylor(:,m_NodDomLoc) = true;
      %Se quita los grados que sean de m�nima restricci�n
      m_DoffNodTaylor = m_DoffNodTaylor&~m_DoffNodMinRest;
      %Matriz de los datos de las condiciones de borde de taylor.
      m_CondBordTaylor = [m_NodDomLoc',m_DoffNodTaylor(:,m_NodDomLoc)'*nroGdl*10.^(ndn-1:-1:0)',...
         zeros(nNodDomLoc,ndn)];
      %Se detecta los nodos que quedan sin restricci�n, quit�ndolos (por ejemplo que haya un nodo con todos
      %sus grados de libertad con m�nima restricci�n).
      m_CondBordTaylor = m_CondBordTaylor(any(m_DoffNodTaylor(:,m_NodDomLoc),1),:);      
      %Se redefine la matriz de datos de condiciones de borde previa anulando los gdl de libertad pisados por
      %la condici�n de borde de Taylor.
      m_DirRestr(m_DoffNodTaylor(:,m_NodRestr)) = 0;
      m_CondBord(:,2) = m_DirRestr'*10.^(ndn-1:-1:0)';
      %Se detecta los grados de libertad que no pisa la condici�n de borde de Taylor.
      m_DoffNodNoTaylor = m_DoffNod&~m_DoffNodTaylor;
      %Se quita los datos de condiciones de borde, los nodos que se quedaron sin restricci�n.      
      m_CondBord = m_CondBord(any(m_DoffNodNoTaylor(:,m_NodRestr),1),:);
   else
      m_CondBordTaylor = [m_NodDomLoc',sum(nroGdl*10.^(0:ndn-1))*ones(nNodDomLoc,1),zeros(nNodDomLoc,ndn)];
   end
   %
   m_CondBord = [m_CondBord;m_CondBordTaylor];

end