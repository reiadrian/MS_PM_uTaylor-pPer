function m_GdlRestMR = f_DetGdlFijoCBMR(m_LinCondRestHom,m_ConecFrontCU,doff,ndoft,ndn)

   ceroDetRel = 1e-3;

   %Matriz con los n�meros de nodos
   m_NodLibres = 1:ndoft/ndn;   
   
   %Se quita los nodos que ya tiene una restricci�n (si un gdl de libertad tiene al menos una
   %restricci�n, se quita ese nodo). Tambi�n se quita los nodos de la frontera de la celda unitaria 
   %(se toma como elementos de frontera los correspondientes a la primera condici�n de 
   %homogenizaci�n)
   m_indNodVal = any(reshape(doff',ndn,[])',2);
   %m_indNodVal(m_ConecFrontCU) = true;
   m_NodFijos = m_NodLibres(m_indNodVal);
   %m_NodLibres = m_NodLibres(~m_indNodVal); 
   %Grados de libertad que pueden ser elegidos
   %m_GdlNodLibres = f_DofElem(m_NodLibres,ndn);
   m_GdlNodFijos = f_DofElem(m_NodFijos,ndn);
   %m_LinCondRestHom = m_LinCondRestHom(:,m_GdlNodLibres);
   %m_LinCondRestHomLibresAbs = abs(m_LinCondRestHom(:,m_GdlNodLibres));
   %Analizar si no vale la pena transformar la matriz m_LinCondRestHomLibresAbs en full.
   m_LinCondRestHomLibresAbs = full(abs(m_LinCondRestHom));
   %A los gdl correspondiente que tienen al menos una restricci�n se le impone NaN para no
   %considerarlos.
   m_LinCondRestHomLibresAbs(:,m_GdlNodFijos) = NaN;
   
   %Determinaci�n de los tres grados de libertad, el nodo al que se impone el gdl con direcci�n 
   %�nica.
   m_GdlRestMR = zeros(3,1);
   %Se basa en tomar el mayor en valor absoluto, y por lo tanto m�s distinto de cero.
   %Es importante que en la primera fila de m_LinCondRestHom est� la ecuaci�n de ux x nx, y en la
   %segunda fila, est� la ecuaci�n correspondiente a los uy x ny, y en que en la �ltima fila est� la
   %ecuaci�n correspondiente a 1/2*(ux * ny + uy * nx).
   [m_Valor,m_IndMaxGdl] = max(m_LinCondRestHomLibresAbs(1:2,:),[],2);
   %En el primero se guarda el gdl en x (es importante este orden por como se define dirgdl).
   %m_GdlRestMR(1) = m_GdlNodLibres(m_IndMaxGdl(1));
   %Se considera que el �ndice de la matriz coincide con el gdl correspondiente.
   m_GdlRestMR(1) = m_IndMaxGdl(1);
   %En el segundo se guarda el gdl en y.
   %m_GdlRestMR(2) = m_GdlNodLibres(m_IndMaxGdl(2));
   m_GdlRestMR(2) = m_IndMaxGdl(2);
   %Se calcula los nodos correspondientes a esos gdl.
   %nod1 = ceil(m_GdlNodLibres(m_IndMaxGdl(1))/ndn);
   %nod2 = ceil(m_GdlNodLibres(m_IndMaxGdl(2))/ndn);
   %Para que eso gdl no se considere m�s ese gdl en la b�squeda del pr�ximo gdl.
   m_LinCondRestHomLibresAbs(3,m_IndMaxGdl) = NaN;
   if m_Valor(1)>=m_Valor(2)
      %C�mo el gdl x result� m�s grande, se adopta este como el gdl �nico, y se busca el tercer gdl
      %con direcci�n y (2). Por eso se pone NaN en los gdl x.     
      dirgdl = 2;
      m_LinCondRestHomLibresAbs(3,1:ndn:ndoft) = NaN;
   else
      %Idem, pero se busca los gdl en x.
      dirgdl = 1;
      m_LinCondRestHomLibresAbs(3,2:ndn:ndoft) = NaN;
   end
   
   %El tercer gdl se busca uno que tenga en la tercer fila, en los grados de libertad dirgdl,
   %el mayor valor, y se verifica que el determinante de 2x2 de los gdl repetidos suficientemente 
   %distinto de cero.
   [~,IndMaxGdl3] = max(m_LinCondRestHomLibresAbs(3,:),[],2);
   %m_GdlRestMR(3) = m_GdlNodLibres(IndMaxGdl3);
   m_GdlRestMR(3) = IndMaxGdl3;
   %Se arma la matriz de condiciones de borde
   m_CBRed = full([m_LinCondRestHom(dirgdl,m_GdlRestMR(dirgdl)),m_LinCondRestHom(dirgdl,m_GdlRestMR(3));...
      m_LinCondRestHom(3,m_GdlRestMR(dirgdl)),m_LinCondRestHom(3,m_GdlRestMR(3))]);
   %Se verifica que el valor relativo del determinante (det(Mat)/norm(Mat(1:2,1))/norm(Mat(1:2,2))
   %sea distinto de cero, usando una tolerancia.
   while (m_CBRed(1)*m_CBRed(4)-m_CBRed(3)*m_CBRed(2))/hypot(m_CBRed(1),m_CBRed(2))/...
         hypot(m_CBRed(3),m_CBRed(4))<ceroDetRel
      m_LinCondRestHomLibresAbs(3,IndMaxGdl3) = NaN;
      [valor3,IndMaxGdl3] = max(m_LinCondRestHomLibresAbs(3,:),[],2);
      %m_GdlRestMR(3) = m_GdlNodLibres(IndMaxGdl3);
      m_GdlRestMR(3) = IndMaxGdl3;
      m_CBRed(2,1) = m_LinCondRestHom(3,m_GdlRestMR(dirgdl));
      m_CBRed(2,2) = m_LinCondRestHom(3,m_GdlRestMR(3));
      if isnan(valor3)
         error(['Determinaci�n de nodos fijos de la CB de M�n. Restr: ',...
            'No se encontr� ning�n nodo con las condiciones buscadas.'])
      end
   end
   fprintf('Grados de libertad elegidos: (Comp x): %d (Comp y): %d (Comp xy): %d\n',m_GdlRestMR) 

end