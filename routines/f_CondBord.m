function [m_LinCond,m_CteCond,m_InvCRR,doff,dofl,doffCondCte] = f_CondBord(e_VG,...
   xx,e_DatSet,m_ConecFront,m_CteMinRestr)

   %Se lee las variables de la estructura e_VG.
   m_CondBordLeid = e_VG.m_CondBord;
   %Ver si no poner m_ConecFront y m_CteMinRestr dentro de las opciones de las condiciones de borde.
   c_OpCondBord = e_VG.c_OpCondBord;
   ndoft = e_VG.ndoft;
   ndn = e_VG.ndn;
   ndime = e_VG.ndime;
   nSet = e_VG.nSet;
   struhyp = e_VG.struhyp;
   
   %Corregir que con la condicion de MR si o si debe ponerse los elementos de frontera, ya que sino el
   %programa funciona sin error y la matriz de CB da NaN sin aviso de donde esta el error.
   %
   sinCond = 0;
   condCte = 1;
   condIgDesp = 2;
   condHomogDef = 3;
   condMedia = 6;
   %Para condiciones de borde variable (segunda frontera)
   condHomogDef2 = 4;
   if iscell(m_ConecFront)
      m_ConecFront2 = m_ConecFront{2};
      m_ConecFront = m_ConecFront{1};
   end
   
   m_CondBanda=[];
   m_CondBord=[];
   i_CondBanda=0;
   i_CondBorde=0;
   CondBand   =0;
   for i_CB= 1:size(m_CondBordLeid,1)
       if m_CondBordLeid(i_CB,2) == 55 || m_CondBordLeid(i_CB,2) == 51|| m_CondBordLeid(i_CB,2) == 15
             i_CondBanda=i_CondBanda+1;
             m_CondBanda(i_CondBanda,:)= m_CondBordLeid(i_CB,:);
             CondBand =1;
       else
             i_CondBorde=i_CondBorde+1;
             m_CondBord(i_CondBorde,:)= m_CondBordLeid(i_CB,:);
       end
   end
   m_CondBordLeid = m_CondBord;
   
   %% Tratamiento de los grados de libertad de las restricciones impuestas
   %Se ordena en la filas ndn y en las columnas la cantidad de nodos con restriccion porque facilita
   %en algunas partes el indexado.
   %Se transforma la indicacion unica de grados de libertad restringidos en columnas separadas.
   m_DirRestr = zeros(ndn,size(m_CondBordLeid,1));
   m_DirRestr(ndn,:) = m_CondBordLeid(:,2);
   for iGdl = 1:ndn-1
      decGdl = 10^(ndn-iGdl);
      m_DirRestrGdl = fix(m_DirRestr(ndn,:)/decGdl);
      m_DirRestr(iGdl,:) = m_DirRestrGdl;
      m_DirRestr(ndn,:) = m_DirRestr(ndn,:)-m_DirRestrGdl*decGdl;
   end   
   %Grados de libertad de los nodos en que se le impuso la restriccion
   m_GdlNodRestr = reshape(f_DofElem(m_CondBordLeid(:,1)',ndn),ndn,[]);

   %Valores que se le pone a los nodos restringidos
   m_ValNodRestr = m_CondBordLeid(:,3:end)';

   %La celdas con las opciones se traspone para coincidir con la forma de ordenar las matrices
   c_OpCondBord = c_OpCondBord';

   %% Grados de libertad globales restringidos (f: fijos) y libres (l)
   if 0
      doff = m_GdlNodRestr(m_DirRestr~=sinCond);
      %Se hace la verificacion de condiciones de borde a nivel de grados de libertad en lugar de
      %nodos para que se pueda ingresar dos lineas de condiciones de borde aplicadas a los mismos
      %nodos, pero grados de libertad distintos. Esto permite que valores aplicados a los grados de
      %libertad sean distintos en caso de periodicidad y a la vez que permite condiciones de
      %periodicidad con los de rigidez.
      if length(doff)~=length(unique(doff))
         error(['Lectura de datos: Condiciones de Borde: Determinacion de matrices: Los grados',...
            ' de libertad con restriccion no deben estar repetidos.'])
      end
      m_esGdlLibre = true(ndoft,1);
      m_esGdlLibre(doff) = false;
      dofl = (1:ndoft)';
      dofl = dofl(m_esGdlLibre);
      nGdlf = length(doff);
   else
      doff = false(ndoft,1);
      %Si tiene algun grado de libertad repetido, no influye en la cantidad de true ingresado en la
      %matriz doff.
      doff(m_GdlNodRestr(m_DirRestr~=sinCond)) = true;
      dofl = ~doff;
      nGdlf = sum(doff);
      %Se hace la verificacion de condiciones de borde a nivel de grados de libertad en lugar de
      %nodos para que se pueda ingresar dos lineas de condiciones de borde aplicadas a los mismos
      %nodos, pero grados de libertad distintos. Esto permite que valores aplicados a los grados de
      %libertad sean distintos en caso de periodicidad y a la vez que permite condiciones de
      %periodicidad con los de rigidez.
      if sum(sum(m_DirRestr~=sinCond))~=nGdlf
         error(['Lectura de datos: Condiciones de Borde: Determinacion de matrices: Los grados',...
            ' de libertad con restriccion no deben estar repetidos.'])
      end
   end
   
   %% Desplazamientos prescriptos impuestos
   %Cuidado hay que forzar que any analice en direccion de las filas, porque si no va funcionar bien en el
   %caso de un problema de un solo grado libertad, ya que en ese caso de un vector columna analiza en
   %direccion de las columnas.
   %#####################################################################################################
   m_IndGdlRestrCondCte = m_DirRestr==condCte;
   [m_RestrLinCondCte,m_RestrCteCondCte,m_GdlRestrCondCte,m_GdlVincCondCte] = f_RestrCondCte(...
      m_GdlNodRestr(:,any(m_IndGdlRestrCondCte,1)),m_ValNodRestr(:,any(m_IndGdlRestrCondCte,1)),...
       m_IndGdlRestrCondCte(:,any(m_IndGdlRestrCondCte,1)));
   %#####################################################################################################
   %Determinacion de que grados de libertad son pertenecientes a la condicion de movimientos rigidos.
   if nargout>4
      %Solo sirve si se utiliza matrices logicas para llevar los grados de libertad libre y fijos.
      %Falta completar para el otro caso.
      doffCondCte = false(ndoft,1);
      doffCondCte(m_GdlNodRestr(m_IndGdlRestrCondCte)) = true;
   end
    
   %% Restricciones de igualdad de desplazamientos
   %#####################################################################################################
   %AA22: agregue e_VG.protype
   m_IndGdlRestrIgDesp = m_DirRestr==condIgDesp;
   [m_RestrLinIgDesp,m_RestrCteIgDesp,m_GdlRestrIgDesp,m_GdlVincIgDesp,m_GdlRestrCteIgDesp] = ...
      f_RestrIgDesp(m_GdlNodRestr(:,any(m_IndGdlRestrIgDesp,1)),...
      m_ValNodRestr(:,any(m_IndGdlRestrIgDesp,1)),m_IndGdlRestrIgDesp(:,any(m_IndGdlRestrIgDesp,1)),...
      xx(:,1:ndime),ndn,nSet,e_DatSet,c_OpCondBord(:,any(m_IndGdlRestrIgDesp,1)),e_VG.protype);
   %#####################################################################################################
   %% Restricciones por homogenizacion de las deformaciones
   m_IndGdlRestrHom = m_DirRestr==condHomogDef;
   if any(any(m_IndGdlRestrHom))
      [m_RestrLinHom,m_RestrCteHom,m_GdlRestrHom,m_GdlVincHom,m_GdlRestrCteHom] = ...
         f_RestrCondHomog(...
         m_GdlNodRestr(:,any(m_IndGdlRestrHom,1)),m_ValNodRestr(:,any(m_IndGdlRestrHom,1)),...
         m_IndGdlRestrHom(:,any(m_IndGdlRestrHom,1)),xx(:,1:ndime),ndime,nSet,ndn,struhyp,e_DatSet,...
         m_ConecFront);
   else
      m_RestrLinHom = [];
      m_RestrCteHom = [];
      m_GdlRestrHom = [];
      m_GdlVincHom = [];
      m_GdlRestrCteHom = [];
   end
   
   %% Restricciones por homogenizacion de las deformaciones
   %m_IndGdlRestrHom = m_DirRestr==condHomogDef2;
   %Corregir lo siguiente para indicar que se impone la segunda condicion de borde de
   %homogenizacion. Tambien corregir el parche siguiente. Hacer que las funciones se le pase un
   %vector con los gdl y los valores impuestos.
   if exist('m_ConecFront2','var')
      %Se impone cualquier valor
      m_IndGdlRestrHom = [true,false;true,true];
   else
      m_IndGdlRestrHom = false(ndn,2);
   end
   if any(any(m_IndGdlRestrHom))
      m_GdlNodRestrHom2 = zeros(ndn,2);
      m_GdlNodRestrHom2(m_IndGdlRestrHom) = 1:3;
      m_ValNodRestr = zeros(ndn,2);
      [m_RestrLinHom2,m_RestrCteHom2,m_GdlRestrHom2,m_GdlVincHom2,m_GdlRestrCteHom2] = ...
         f_RestrCondHomog(...
         m_GdlNodRestrHom2,m_ValNodRestr,...
         m_IndGdlRestrHom,xx(:,1:ndime),ndime,nSet,ndn,struhyp,e_DatSet,m_ConecFront2,...
         m_CteMinRestr);
      %Se hace una iteracion por lo menos para asegurar que la parte correspondiente a la segunda
      %condicion de borde de homogenizacion esta bien condicionada (determinante distinto de cero
      %alejado de cero). Igual puede fallar al combinarse con otras condiciones de borde.
      %Para no tener que hacer este ensamblaje, habria que elegir gdl arbitrarios para cada fila de
      %la matriz (la posicion de las filas no cambian el resultado), y asi poder ensamblar sin
      %conocer los gdl fijos de esta condicion, y luego tomar esas filas para determinar la 
      %matriz m_LinCondRestHom2.
      m_LinCondRestHom2 = sparse(m_GdlRestrHom2,m_GdlVincHom2,m_RestrLinHom2,3,ndoft);
      %Se asume que la conectividad de la frontera indicada en la primera es la frontera de la
      %malla.
      m_GdlRestMR = f_DetGdlFijoCBMR(m_LinCondRestHom2,m_ConecFront2,doff,ndoft,ndn);    %Aca estaba m_ConecFront en lugar de m_ConecFront2, pero no se usa para nada asi que no influye
      doff(m_GdlRestMR) = true;
      dofl(m_GdlRestMR) = false;
      nGdlf = nGdlf+3;
      m_GdlRestrHom2 = repmat(m_GdlRestMR,length(m_GdlVincHom2)/3,1);
      %condMat = cond(full(m_LinCondRestHom2(:,m_GdlNodRestr(m_IndGdlRestrHom))));
      condMat = condest(m_LinCondRestHom2(:,m_GdlRestMR));
      fprintf('Numero de condicion %f\n',condMat);
      %while 1/condMat<1e8*eps(1)    %condMat>1e-8/eps(1)
      %end
   else
      m_RestrLinHom2 = [];
      m_RestrCteHom2 = [];
      m_GdlRestrHom2 = [];
      m_GdlVincHom2 = [];
      m_GdlRestrCteHom2 = [];
   end
   
%    [m_RestrLinCondCte,m_RestrCteCondCte,m_GdlRestrCondCte,m_GdlVincCondCte] = f_RestrBanda(...
%       m_GdlNodRestr(:,any(m_IndGdlRestrCondCte)),m_ValNodRestr(:,any(m_IndGdlRestrCondCte)),...
%        m_IndGdlRestrCondCte(:,any(m_IndGdlRestrCondCte)));
   
  
   %% restriccion de condicion de borde de la media del campo de desplazmiento
   m_IndGdlRestrMedia = m_DirRestr==condMedia;
   %#####################################################################################################
   [m_RestrLinCondMedia,m_RestrCteCondMedia,m_GdlRestrCondMedia,m_GdlVincCondMedia,m_GdlRestrCteCondMedia] = ...
      f_RestrMedia(m_GdlNodRestr(:,any(m_IndGdlRestrMedia,1)),...
      m_ValNodRestr(:,any(m_IndGdlRestrMedia,1)),m_IndGdlRestrMedia(:,any(m_IndGdlRestrMedia,1)),...
      ndn,nSet,e_DatSet);
    %#####################################################################################################
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   %% Matrices de restriccion
   m_LinCond = sparse([m_GdlRestrCondCte;m_GdlRestrIgDesp;m_GdlRestrHom;m_GdlRestrHom2;m_GdlRestrCondMedia],...
      [m_GdlVincCondCte;m_GdlVincIgDesp;m_GdlVincHom;m_GdlVincHom2;m_GdlVincCondMedia],...
      [m_RestrLinCondCte;m_RestrLinIgDesp;m_RestrLinHom;m_RestrLinHom2;m_RestrLinCondMedia],...
      ndoft,ndoft);
   m_CteCond = sparse([m_GdlRestrCondCte;m_GdlRestrCteIgDesp;m_GdlRestrCteHom;m_GdlRestrCteHom2;...
      m_GdlRestrCteCondMedia],...
      1,[m_RestrCteCondCte;m_RestrCteIgDesp;m_RestrCteHom;m_RestrCteHom2;m_RestrCteCondMedia],ndoft,1);
 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   if CondBand
      
      %% Desplazamientos prescriptos impuestos por bandas de localizacion
      %
      %              Nodo 1 -------  Nodo 2       Matriz Boundary Conditions
      %                 |              |
      %                 |     Banda j  |              j    55   Nodo1  Nodo2
      %                 |              |              j    55   Nodo3  Nodo4
      %              Nodo 3 -------  Nodo 4           j    55   Nodo5  Nodo6
      %                 |              |              j    55   Nodo7  Nodo8
      %                 |    Banda j   |
      %                 |              |
      %              Nodo 5 -------  Nodo 6
      %                 |              |
      %                 |    Banda j   |
      %                 |              |
      %              Nodo 7 -------  Nodo 8
      %
      
      nband = max(m_CondBanda(:,1));
      CondBand=struct('matriz',[]);
      for i_band=1:nband
          CondBand(i_band).matriz = m_CondBanda(find(m_CondBanda(:,1)==i_band),:);
          % CondBand{i_band}.matrix= m_CondBanda(find(m_CondBanda(:,1)==i_band),:);
      end
      
      for i_band=1:nband
          node1= CondBand(i_band).matriz(1,3);
          node2= CondBand(i_band).matriz(1,4);
          doff1= [((node1-1)*ndn+1) ; node1*ndn];
          doff2= [((node2-1)*ndn+1) ; node2*ndn];
          nCond_Band= size(CondBand(i_band).matriz,1);
          
          for i_eq= 2:nCond_Band;
              
              node3= CondBand(i_band).matriz(i_eq,3);
              node4= CondBand(i_band).matriz(i_eq,4);
              doff3= [((node3-1)*ndn+1) ; node3*ndn];
              doff4= [((node4-1)*ndn+1) ; node4*ndn];
              
              doff(doff3)=1;
              nGdlf=nGdlf+2;
              m_LinCond(doff3(1),doff3(1))=1;m_LinCond(doff3(1),doff4(1))=-1;...
                  m_LinCond(doff3(1),doff2(1))=1;m_LinCond(doff3(1),doff1(1))=-1;
              
              m_LinCond(doff3(2),doff3(2))=1;m_LinCond(doff3(2),doff4(2))=-1;...
                  m_LinCond(doff3(2),doff2(2))=1;m_LinCond(doff3(2),doff1(2))=-1;
          end
      end
      dofl = ~doff;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % Condensacion para que en funcion de los grados de libertad sin restriccion
   %m_InvCRR = inv(m_LinCond(doff,doff));
   m_InvCRR = m_LinCond(doff,doff)\speye(nGdlf,nGdlf);
   %Ver si es mas lento asi para matrices sparse que hacer dos veces la inversa.
   %m_LinCond = m_LinCond(doff,doff)\[-m_LinCond(doff,dofl),m_CteCond(doff)];
   %m_LinCond = m_LinCond(doff,doff)\-m_LinCond(doff,dofl);
   m_LinCond = m_InvCRR*-m_LinCond(doff,dofl);
   %Vector restricciones constantes de los grados de libertad fijos (doff).
   %m_CteCond = m_LinCond(:,end);
   %m_LinCond = m_LinCond(:,1:end-1);
   
   %Se incorpora esta opcion porque se observo en el caso que las matrices de condiciones de borde que son
   %llenas, como ejemplo en el caso de la media de temperatura (la temperatura de todos los nodos dependen
   %entre si), hace que la matriz de rigidez global condensada tenga las caracteristicas de una matriz llena.
   %Esto ocurre por el esquema de Condensacion directa utilizada para la resolucion del sistema:
   %KT = KT(dofl,dofl)+m_LinCond'*KT(doff,dofl)+KT(dofl,doff)*m_LinCond+m_LinCond'*KT(doff,doff)*m_LinCond;
   %donde principalmente aporta a la perdida de la sparsidad el termino m_LinCond'*KT(doff,doff)*m_LinCond.
   %Que la matriz de rigidez global condensada tenga caracteristica de llena y se almacene en una esquema de
   %matrices sparse hace que funcione mas lento (aproximadamente la mitad) y ocupe mas espacio
   %(aproximadamente el doble, siendo esto lo peor porque puede hacer que swapee, es decir que use la memoria
   %virtual y en consecuencia el disco), que en el caso de utilizar directamente matrices llenas.
   %Al hacer las matrices de condiciones de borde full, hace que la matriz rigidez global condensada sea full,
   %aunque la matriz de rigidez completa (no condensada) sea sparse.
   %Por ahora se convierte la matrices a full, habria que ver si no organizarlo de otra forma para ensamblarlo
   %directamente como matrices full.
   if e_VG.matCBFull
      m_LinCond = full(m_LinCond);
      m_CteCond = full(m_CteCond);
   end

end
 
%%
function esInt = f_esPtoIntElem2D(m_CoordNod,PtoInt,eltype)
   
   %Se asume que los nodos de las conectividades estan ordenadas en sentido antihorario, es decir la
   %lista de coordenadas en m_CoordNod estan ordenadas de esa forma. Esto es importante para la
   %verificacion que el punto sea interior.
   %Solo sirve cuando el elemento tiene forma convexa en la malla (en coordenadas globales).
   if eltype==16
       m_CoordNodA = zeros(size(m_CoordNod));
       m_CoordNodA(1:2:size(m_CoordNod,1),:) = m_CoordNod(1:4,:);
       m_CoordNodA(2:2:size(m_CoordNod,1),:) = m_CoordNod(5:8,:);
       m_CoordNod = m_CoordNodA;
   end
   m_DifCoord = [diff(m_CoordNod);m_CoordNod(1,:)-m_CoordNod(end,:)];
   m_DifPtoInt = bsxfun(@minus,PtoInt,m_CoordNod);
   esInt = all((m_DifPtoInt(:,1).*m_DifCoord(:,2)-m_DifPtoInt(:,2).*m_DifCoord(:,1))<=0);
   
end

%%
function [m_RestrLin,m_RestrCte,m_GdlRestr,m_GdlVinc,m_GdlRestrCte] = f_RestrIgDesp(...
   m_GdlNodRestr,m_ValNodRestr,m_IndGdlRestr,xx,ndn,nSet,e_DatSet,c_OpCondBord,protype)

   % Restricciones de igualdad de desplazamientos
   %En forma general no se conoce la cantidad de grados de libertad que van estar vincunlados con
   %la igualdad de desplazamiento, ya que no se conoce el tipo de elemento donde coincide la
   %posicion que se impone la misma.
   m_GdlVinc = [];
   m_GdlRestr = [];
   m_RestrLin = [];
   %
   if ~isempty(m_IndGdlRestr)
      nNodRestr = size(m_GdlNodRestr,2);   
      for iNodRestr = 1:nNodRestr
         ptoEncon = false;
         %Se asume que si se agrega datos en las opciones es porque las coordenadas a la que esta vinculado el
         %nodo se ingresa por ahi.
         %Se asume que c_OpCondBord tiene un valor por nodo asignado.
         s_OpCondBord = c_OpCondBord{:,iNodRestr};
         nousarOp = isempty(s_OpCondBord);
         if nousarOp
            m_Pos = m_ValNodRestr(:,iNodRestr);
         else
            m_Pos = textscan(c_OpCondBord{:,iNodRestr},'(%f,%f)','CollectOutPut',true);
            m_Pos = m_Pos{1}';
         end
         %##################################################################
         if protype==1 %AA22: Agregue 2022 
             gdl_periodico=sum(double(m_IndGdlRestr(:,iNodRestr)));
             if gdl_periodico == 2 || gdl_periodico==3
                 m_Pos = m_Pos(1:2);
             elseif gdl_periodico == 1
                 m_Pos = m_Pos(3:4);
             end
         end %AA22: Agregue 2022 
         %##################################################################
         for iSet =  1:nSet
            eltype = e_DatSet(iSet).e_DatElem.eltype;
            conec = e_DatSet(iSet).conec;
            m_DofElem = e_DatSet(iSet).m_DofElem;
            %Cantidad de elementos en el Set.
            nElem = e_DatSet(iSet).nElem;
            for iElem = 1:nElem
               m_CoordNodElem = f_CoordElem(xx,conec(iElem,:))';
               if f_esPtoIntElem2D(m_CoordNodElem,m_Pos',eltype)
                  switch eltype
                     case {2,32}
                        %Triangulo de tres nodos con funciones de forma (para aprox. de coord.) lineales
                        m_FunFormX = [m_CoordNodElem,ones(3,1)]'\[m_Pos;1];
                     case {4,8,20,31,108}
                        %Cuadrangulo de cuatro nodos con funciones forma (para aprox. de coord.) 
                        %bilineales
                        m_FunFormX = [m_CoordNodElem,m_CoordNodElem(:,1).*m_CoordNodElem(:,2),...
                           ones(4,1)]'\[m_Pos;m_Pos(1)*m_Pos(2);1];
                     case {16} %AA22: Agregue 2022
                        %Cuadrangulo de ocho nodos con funciones forma (para aprox. de coord.) 
                        %bicuadraticas
                        A=[m_CoordNodElem,m_CoordNodElem(:,1).*m_CoordNodElem(:,2),...
                           ones(8,1)]';
                        B=[m_Pos;m_Pos(1)*m_Pos(2);1];
                        m_FunFormX = all(A(:,1:8)==B)';
%                         m_FunFormX = [m_CoordNodElem,m_CoordNodElem(:,1).*m_CoordNodElem(:,2),...
%                            ones(8,1)]'\[m_Pos;m_Pos(1)*m_Pos(2);1];
                                      %AA22: Agregue 2022
                     otherwise
                        error(['Lectura de datos: Condiciones de borde: Determinacion de matrices: ',...
                           'Restricciones de igualdad de desplazamiento: Tipo de elemento no ',...
                           'implementado.']);
                  end
                  %Se guarda como vector para usar el ensamblaje de la funcion sparse.
                  %Se considera que los grados de libertad del elemento estan aproximados con las mismas
                  %funciones de forma.
                  m_GdliNodRestr = m_GdlNodRestr(m_IndGdlRestr(:,iNodRestr),iNodRestr);
                  nGdlNodRestr = length(m_GdliNodRestr);
                  %m_GdlElem = reshape(f_DofElem(conec(iElem,:),ndn),ndn,[]);
                  m_GdlElem = reshape(m_DofElem(:,iElem),ndn,[]);
                  %
                  m_GdlRestr = [m_GdlRestr;repmat(m_GdliNodRestr,size(conec,2)+1,1)];
                  m_GdlVinc = [m_GdlVinc;m_GdliNodRestr;reshape(m_GdlElem(m_IndGdlRestr(:,iNodRestr),...
                     :),[],1)];
                  m_RestrLin = [m_RestrLin;ones(nGdlNodRestr,1);reshape(repmat(-m_FunFormX,1,...
                     nGdlNodRestr)',[],1)];
                  %Se considera que el desplazamiento se evalua silo con la funcion de forma de un
                  %elemento. Cuando coincide la posicion con un lado contiguo de dos elementos se evalua
                  %solo sobre uno de ellos. (Esto no es correcto con mallas compuestas?).
                  %Varias restricciones pueden estar sobre los mismos elementos.
                  ptoEncon = true;
                  break
               end
            end
         end
         if ~ptoEncon
            warning(['Lectura de datos: Condiciones de borde: Determinacion de matrices: ',...
               'Restricciones de igualdad de desplazamiento: El punto (x,y)=(%g,%g) vinculado ',...
               'a un grado de libertad restringido no fue encontrado sobre la malla.'],m_Pos)
         end
      end
      %Este no se utiliza pero deberia devolver los valores que se ensamblan en la matriz de constante.
      %Es para si se quiere que los desplazamientos vinculados tenga una diferencia en su valor. En el caso
      %que las coordenadas se ingresan como opciones el dato ingresado en el valor de desplazamiento
      %prescripto representa la diferencia entre los desplazamientos que se impone.
      m_GdlRestrCte = m_GdlNodRestr(m_IndGdlRestr);
      %Para asegurar que con un solo grado de libertad sea un vector columna.
      m_GdlRestrCte = m_GdlRestrCte(:);
      if nousarOp
         m_RestrCte = zeros(length(m_GdlRestrCte),1);
      else         
         m_RestrCte = m_ValNodRestr(m_IndGdlRestr);
         m_RestrCte = m_RestrCte(:);
      end
   else
      %Para evitar problemas en la indexacion  
      m_GdlRestrCte = [];
      m_RestrCte = [];
   end
   
end

%%
function [m_RestrLin,m_RestrCte,m_GdlRestr,m_GdlVinc] = f_RestrCondCte(m_GdlNodRestr,...
   m_ValNodRestr,m_IndGdlRestr)
   
   % Restricciones debidas a desplazamientos prescriptos
   %Para que funcione correctamente la concatenacion cuando las matrices son vAcaos pero con dimension, se
   %fuerza que sean vAcaas y sin dimension. Ver si se puede hacer automaticamente, es decir que salga de la
   %indexacion (cuidado que no es tan simple y recordar que deberia funcionar para distintos Numeros de
   %grados de libertad).
   %###########################################################################
   %Caso modelo sin restriccion (MODEL=14) y con metodo de la penalidad (MODEL=15)
   %Al imponer la periodicidad de las poro presiones necesito el siguiente formato de archivo de entrada
   %-Nodo Restriccion_ux Restriccion_uy Restriccion_p Valor_ux Valor_uy Valor_p1 Valor_p2-
   %Donde 
   %si Restriccion_ux/Restriccion_uy = 0 --> No existe restriccion (libre)
   %si Restriccion_ux/Restriccion_uy = 1 --> Valor_ux/Valor_uy=valor de desplazamiento finito 
   %si Restriccion_ux/Restriccion_uy = 2 --> Valor_ux=coord. x del nodo 'opuesto' /Valor_uy=coord. y del nodo 'opuesto'
   %En cambio para p deberia funcionar asi
   %Restriccion_p = 0 --> No existe restriccion (libre)
   %Restriccion_p = 1 --> Valor_p1 = valor de poro presion finito. BASTA
   %CON UN SOLO VALOR, EL Valor_p2 no haria falta y es redundante.
   %Restriccion_p = 1 --> Valor_p1 =coord. x del nodo 'opuesto'/ Valor_p2=coord. y del nodo 'opuesto'
   %En esta funcion que existan 4 valores (Valor_ux Valor_uy Valor_p1 Valor_p2)
   %conduce a que 'm_ValNodRestr'(4filas,ncol) sea de dimension mayor a 'm_IndGdlRestr' (3filas,ncol)
   %y esto introduce error en m_RestrCte (y por ende en todo el
   %procedimiento de calculo posterior)
   %De esta manera los 2 modelos solo presentan el
   %problema en esta funcion de que m_ValNodRestr(4filas,ncol)>m_IndGdlRestr(3filas,ncol) en numero
   %de filas. Para resolver el problema se adopta que: size(m_ValNodRestr,1) > 3
   %m_ValNodRestr=m_ValNodRestr(1:end-1,:) que RESOLVERIA el problema. QUEDA
   %POR COMPROBAR
   
   %Caso modelo con metodo de multiplicadores de Lagrange (ML) (MODEL=16)
   %Al imponer la periodicidad de las poro presiones necesito el siguiente formato de archivo de entrada
   %-Nodo Restriccion_ux Restriccion_uy Restriccion_p Restriccion_ML Valor_ux Valor_uy Valor_p1 Valor_p2-
   %Donde 
   %si Restriccion_ux/Restriccion_uy = 0 --> No existe restriccion (libre)
   %si Restriccion_ux/Restriccion_uy = 1 --> Valor_ux/Valor_uy=valor de desplazamiento finito 
   %si Restriccion_ux/Restriccion_uy = 2 --> Valor_ux=coord. x del nodo 'opuesto' /Valor_uy=coord. y del nodo 'opuesto'
   %En cambio para p deberia funcionar asi
   %Restriccion_p = 0 --> No existe restriccion (libre)
   %Restriccion_p = 1 --> Valor_p1 = valor de poro presion finito. BASTA
   %CON UN SOLO VALOR, EL Valor_p2 no haria falta y es redundante.
   %Restriccion_p = 2 --> Valor_p1 =coord. x del nodo 'opuesto'/ Valor_p2=coord. y del nodo 'opuesto'
   %En esta funcion que existan 4 valores (Valor_ux Valor_uy Valor_p1 Valor_p2)
   %conduce a que 'm_ValNodRestr'(4filas,ncol) sea de dimension mayor a 'm_IndGdlRestr' (3filas,ncol)
   %y esto introduce error en m_RestrCte (y por ende en todo el
   %procedimiento de calculo posterior)
   %Finalment para la restriccion del multiplicador entiendo siempre debe
   %estar LIBRE. Entonces 
   %Restriccion_ML = 0 --> No existe restriccion (libre)
   %Por esta razon no hace falta un Valor_M. SIEMPRE DEBERIA TOMAR Restriccion_ML = 0
   %De esta manera para este modelo (dado lo ultimo) solo presenta el
   %problema en esta funcion de que m_ValNodRestr(5filas,ncol)=m_IndGdlRestr(4filas,ncol) en numero
   %de filas pero solo interasan los valores hasta la fila 4!
   %Para resolver el problema se CONSIDERA que: size(m_ValNodRestr,1) > 3
   %m_ValNodRestr=m_ValNodRestr(1:end-1,:) que RESOLVERIA el problema. QUEDA POR COMPROBAR
   %Para los modelos 14,15 y 16 en la MICRO-ESCALA size(m_ValNodRestr,1) > 3
   %Para los modelos MACROESCALA la condicion periodica no interesa
   %entonces las filas de m_ValNodRestr = 3 SIEMPRE
   %###########################################################################
   if size(m_ValNodRestr,1) > 3
       m_ValNodRestr=m_ValNodRestr(1:end-1,:);%QUEDA POR COMPROBAR
   end
   %###########################################################################
   if ~isempty(m_IndGdlRestr)
      m_GdlRestr = m_GdlNodRestr(m_IndGdlRestr);
      m_GdlVinc = m_GdlRestr;
      m_RestrCte = m_ValNodRestr(m_IndGdlRestr);
      %Para asegurar que sean columnas estas matrices para cualquier cantidad de grados de libertad por nodo
      %(cuando ndn=1, las matrices que devuelve esta indexacion es filas).
      %Esto podria evitar el problema de matrices vAcaas sin las dimensiones correctas, ver que el problema es
      %el mismo.
      m_GdlRestr = m_GdlRestr(:);
      m_GdlVinc = m_GdlVinc(:);
      m_RestrCte = m_RestrCte(:);
      %
      m_RestrLin = ones(length(m_RestrCte),1);
   else
      m_GdlRestr = [];
      m_GdlVinc = [];
      m_RestrCte = [];
      m_RestrLin = [];      
   end

end

%%
function [m_C,m_RestrCte,m_GdlRestr,m_GdlVinc,m_GdlRestrCte] = f_RestrCondHomog(...
   m_GdlNodRestr,m_ValNodRestr,m_IndGdlRestr,xx,ndime,nSet,ndn,struhyp,e_DatSet,m_ConecFront,m_CteMinRestr)

    %           4   5
    %  11 o-----o---o
    %     |         |
    %     |         |
    %   6 o         |
    %     |         o 15
    %     |         |
    %   1 o--o------o 9
    %        8
    %
    % conec : lista de conectividades para pseudo-elementos de barra en la frontera, 
    % orden en sentido anti-horario. Ejemplo:
    % conec = [ 1 8 ; 8 9 ; 9 15 ; 15 5 ; 5 4 ; 4 11 ; 11 6 ; 6 1 ]
    %
    % xx    : lista de coordenadas de todos los nodos de la malla
    % xx    = [ x1 y1 z1 ; x2 y2 z2 ; x3 y3 z3 ; x4 y4 z4 ; ... ]

   %m_ValNodRestr: no se utiliza en esta restriccion lineal.
   %Se tiene 3 ecuaciones de independientes de la restriccion por homogenizacion de la deformacion
   %en un problema de dos dimensiones
   %ngdLDep = 3;
   %Se agrega la condicion de borde antisimetrica 1/2*(ux*ny-uy*nx)=0 a las condiciones sobre toda la
   %frontera.
   %En el caso de ser homogeneas imponer (si no es, Tambien seria necesario imponer 1/2*(uy*nx-ux*ny) =
   %-1/2*(ux*ny-uy*nx)), la condicion 1/2*(ux*ny-uy*nx)=0 junta a las otras es equivalente a imponer u * v no
   %simetrico.
   ngdLDepMax = 4;
   
   %% Determinacion de los grados de libertad dependientes (impuestos como restriccion)
   m_GdlDep = m_GdlNodRestr(m_IndGdlRestr);
   %Para asegurar que m_GdlDep sea siempre un vector columna, cualquiera sea la cantidad de Numero de grados
   %de libertad por nodo (cuando tiene un solo grado de libertad, m_IndGdlRestr es un vector fila, y por lo
   %tanto al indexar en m_GdlNodRestr, que es un vector fila devuelve, devuelve un vector fila).
   m_GdlDep = m_GdlDep(:);
   ngdLDep = length(m_GdlDep);
   %ngdLDep = sum(sum(m_IndGdlRestr));
  
   %% Comprobacion
   %m_ElType = [e_DatSet.eltype]';
   m_ElType = arrayfun(@(x)x.e_DatElem.eltype,e_DatSet);
   if any(~(m_ElType==2|m_ElType==4|m_ElType==8|m_ElType==20|m_ElType==31|m_ElType==32|m_ElType==108))
      error(['Lectura de datos: Condiciones de borde: Determinacion de matrices: Restricciones ',...
         'por homogenizacion de la deformacion: Tipo de elemento no implementado.'])
   end
   if ndime~=2
      error(['Lectura de datos: Condiciones de borde: Determinacion de matrices: Restricciones ',...
         'por homogenizacion de la deformacion: Implementado silo para problemas 2D.'])
   end
   if ngdLDep>ngdLDepMax
      error(['Lectura de datos: Condiciones de borde: Determinacion de matrices: Restricciones ',...
         'por homogenizacion de la deformacion: esta implementado para que tenga a lo sumo cuatro ',...
         'restricciones de este tipo.'])
   end
   
   %% Calculo de la normal al elemento
   %(esto es valido para elementos de frontera lineales, con n=cte a lo largo de la longitud de 
   %cada elemento)
   %Vector Tangente (tiene la longitud del elemento)
   tvector = xx(m_ConecFront(:,2),:)-xx(m_ConecFront(:,1),:);
   %Longitud del elemento lineal
   le = sqrt(sum(tvector.^2,2));
   %Vector Tangente normalizado
   tvector = bsxfun(@rdivide,tvector,le);
   %Vector Normal normalizado
   %Aca se define que los elementos de frontera se ingresa en sentido horario para que la normal
   %quede hacia afuera pero como las restricciones de homogenizacion son ecuaciones homogeneas no
   %importa el signo que tenga las componentes de la normal ya que da el mismo resultado (es decir
   %no interesa en que direccion se ingrese los elementos mientras que sea la misma en todos ellos).
   nvector = [tvector(:,2),-tvector(:,1)];

   %% Determinacion de la matriz C global de homogenizacion
   [nelFront,nNodElFront] = size(m_ConecFront);   

   %La cantidad de columnas esta dictado por la cantidad grados de libertad del elemento, por la cantidad de
   %grados de libertad de cada nodo, correspondiente a las funciones de forma de interpolacion de la variable
   %del problema.
   m_C = zeros(ngdLDep,nNodElFront*ndn,nelFront);
   m_GdlVinc = zeros(ngdLDep,nNodElFront*ndn,nelFront);
   
   %Se adopta en forma arbitraria que el primer grado de libertad donde se le impone restriccion de
   %este en el programa tiene la ecuacion de homogenizacion de las componentes xx del tensor, la
   %segunda las componentes yy, y el tercero tiene las componentes xy.
   %Igual el orden no influye en que el sistema esta mal condicionada (cercano a matriz singular).
   m_GdlRestr = repmat(m_GdlDep,[1,nNodElFront*ndn,nelFront]);
         
   %La ubicacion de los valores en la matriz m_C de los grados de libertad depende de como se defina las
   %matrices de funciones de forma del campo que se esta interpolando. Para tener en cuento esto se lo mide
   %con la hipotesis estructural.
   switch struhyp 
      case {1,2,20} %Small deformation Plane strain and Plane stress, large deformations with plane deformation.
         m_IndNPosFFRes1 = 1:ndn:nNodElFront*ndn;
         m_IndNPosFFRes2 = 2:ndn:nNodElFront*ndn;
      case 30   %Thermal problems
         m_IndNPosFFRes1 = 1:ndn:nNodElFront*ndn;
         m_IndNPosFFRes2 = 1:ndn:nNodElFront*ndn;
      otherwise
         error('Lectura de datos: Condiciones de Borde: minima restriccion: Estado no disponible.'); 
   end
   %
   for iElemFront = 1:nelFront
    
      m_DofElemFront = f_DofElem(m_ConecFront(iElemFront,:),ndn);
    
      m_GdlVinc(:,:,iElemFront) = repmat(m_DofElemFront',ngdLDep,1);
    
      %Integral de la funcion de forma sobre el borde del elemento (sobre el elemento de barra).
      %En el caso de elementos lineales la integral es la mitad de la longitud del elemento
      %(Triangulo de base le y altura 1).
      intN = le(iElemFront)/2;
      
      %Para que funcione si se elije imponer menos ecuaciones de esta condicion.
      %Ver que el orden cual restriccion corresponde a Numero de grado de libertad elegido en el archivo de
      %datos es medio arbitrario ya que depende del orden en el que viene en la matriz de m_GdlDep.
      %Lo que si es determinado es que se impone una condicion de este tipo, se impone ux*nx=0, si se
      %impone dos condiciones, se impone ux*nx=0 y uy*ny=0. Si se impone tres condiciones, ademas de las dos
      %anteriores se impone ux*ny+uy*nx. Si hay cuatro condiciones impone las otras tres mas ux*ny-uy*nx.
      if ngdLDep>0
         %ux*nx=0
         m_C(1,m_IndNPosFFRes1,iElemFront) = nvector(iElemFront,1)*intN;
      end
      if ngdLDep>1
         %uy*ny=0
         m_C(2,m_IndNPosFFRes2,iElemFront) = nvector(iElemFront,2)*intN;
      end
      %%condicion simetrica.
      %Aca se elimina el 1/2 asumiendo que es homogonea la condicion de borde.
      if ngdLDep>2
         %ux*ny+uy*nx
         m_C(3,m_IndNPosFFRes1,iElemFront) = nvector(iElemFront,2)*intN;
         m_C(3,m_IndNPosFFRes2,iElemFront) = nvector(iElemFront,1)*intN;
      end
      %%condicion antisimetrica.
      %Aca se elimina el 1/2 asumiendo que es homogonea la condicion de borde.
      if ngdLDep>3
         %ux*ny-uy*nx
         m_C(4,m_IndNPosFFRes1,iElemFront) = nvector(iElemFront,2)*intN;
         m_C(4,m_IndNPosFFRes2,iElemFront) = -nvector(iElemFront,1)*intN;
      end
         
   end
   
   %Se vectoriza las matrices para usarla con la funcion sparse.
   m_GdlRestr = m_GdlRestr(:);
   m_GdlVinc = m_GdlVinc(:);
   m_C = m_C(:);
   
   %Este no se utiliza pero deberia devolver los valores que se ensamblan en la matriz de constante.
   m_GdlRestrCte = m_GdlDep;
   if exist('m_CteMinRestr','var')
      %m_CteMinRestr: Vector conteniendo las ctes. de minima restriccion ordenada segun las componentes xx, yy, y
      %xy.
      m_RestrCte = m_CteMinRestr(:);
   else
      m_RestrCte = zeros(length(m_GdlRestrCte),1);
   end

end

%%
function [m_RestrLin,m_RestrCte,m_GdlRestr,m_GdlVinc,m_GdlRestrCte] = f_RestrMedia(m_GdlNodRestr,...
   m_ValNodRestr,m_IndGdlRestr,ndn,nSet,e_DatSet)

   %nGdLDep = 2;
   %La fila que se elija para insertar en la matriz m_LinCond la ecuaciones de las restricciones, no tienen
   %importancia, mientras que no se repita la fila en dos ecuaciones distintas. Por ello se puede elegir
   %cualquier grado de libertad restringido, ya sea en x o y, para la ecuacion de restriccion de la media 
   %sobre la componente x o y del desplazamiento.
   %Considerando esto se podria hacer m_GdlDep = m_GdlNodRestr(m_IndGdlRestr);, pero tiene el inconveniente
   %que dependiendo como fueron elegidos los grados de libertad dependientes (orden), se le puede asignar 
   %distinto segun el caso (por ejemplo, la componente x se asigna la ecuacion de x o la de y segun
   %como es el orden de los grados de libertad).
   %Para evitar, se asume que vienen dos grados de libertad fijados, en cualquier orden, pero se
   %inserta primero el grado de libertad en x, que va definir la posicion de la ecuacion x, y luego el grado
   %de libertad y que define la fila de la ecuacion y.
   %###########################################################################################
   if ndn==4 %AA: add
   m_GdlDep = [m_GdlNodRestr(1,m_IndGdlRestr(1,:));m_GdlNodRestr(2,m_IndGdlRestr(2,:));...
               m_GdlNodRestr(3,m_IndGdlRestr(3,:));m_GdlNodRestr(4,m_IndGdlRestr(4,:))]; %AA: add
   elseif ndn==3 %AA: add
   m_GdlDep = [m_GdlNodRestr(1,m_IndGdlRestr(1,:));m_GdlNodRestr(2,m_IndGdlRestr(2,:));...
               m_GdlNodRestr(3,m_IndGdlRestr(3,:))]; %AA: add
   elseif ndn==2
   m_GdlDep = [m_GdlNodRestr(1,m_IndGdlRestr(1,:));m_GdlNodRestr(2,m_IndGdlRestr(2,:))];
   elseif ndn==1
      m_GdlDep = m_GdlNodRestr(1,m_IndGdlRestr(1,:));
   else
      error(['Lectura de datos: Condiciones de borde: Determinacion de matrices: Restricciones ',...
         'de media: No esta implementado para que tenga solo dos restricciones de este tipo.'])
   end
   %###########################################################################################
   nGdLDep = length(m_GdlDep);
   if ~isempty(m_GdlDep)
%       if length(m_GdlDep)~=nGdLDep
%          error(['Lectura de datos: Condiciones de borde: Determinacion de matrices: Restricciones ',...
%             'de media: esta implementado para que tenga solo dos restricciones de este tipo.'])
%       end
      %
      c_RestrLin = cell(nSet,1);
      c_Fil = cell(nSet,1);
      c_Col = cell(nSet,1);

      for iSet = 1:nSet
         %
         nElem = e_DatSet(iSet).nElem;
         dofpe = e_DatSet(iSet).e_DatElem.dofpe;
         eltype = e_DatSet(iSet).e_DatElem.eltype;
         wg = e_DatSet(iSet).e_DatElem.wg;
         npg = e_DatSet(iSet).e_DatElem.npg;
         m_DofElem = e_DatSet(iSet).m_DofElem;
         m_FF = e_DatSet(iSet).m_FF;
         m_DetJT = e_DatSet(iSet).m_DetJT;
         %
         m_Col = reshape(repmat(m_DofElem(:)',nGdLDep,1),[],1);
         m_Fil = reshape(repmat(m_GdlDep,dofpe*nElem,1),[],1);
         %
         m_N = zeros(nGdLDep,dofpe,nElem);
         switch eltype
            case {2,4,8,32,108}
               for iElem = 1:nElem
                  m_pesoPG = m_DetJT(:,iElem).*wg;
                  m_Ne = zeros(nGdLDep,dofpe);
                  for iPG = 1:npg
                     %m_Ne = m_Ne+m_FF(:,:,iPG,iElem)*m_pesoPG(iPG);
                     m_Ne = m_Ne+m_FF(:,:,iPG)*m_pesoPG(iPG);
                  end
                  m_N(:,:,iElem) = m_Ne;
               end
            otherwise
               error(['Lectura de datos: Condiciones de borde: Determinacion de matrices: Restricciones ',...
                  'de media: Elemento Finito no definido.'])
         end
         %
         c_RestrLin{iSet} = m_N(:);
         c_Col{iSet} = m_Col;
         c_Fil{iSet} = m_Fil;      
      end
      m_RestrLin = cat(1,c_RestrLin{:});
      m_GdlVinc = cat(1,c_Col{:});
      m_GdlRestr = cat(1,c_Fil{:});

      %Se fija que el termino de la restriccion sea cero, ignorando la que viene del archivo de datos (ver si
      %imponer de la que viene del archivo de datos).
      %m_RestrCte = zeros(m_GdlDep,1);
      %Se asume (ver arriba) que, primero, se impone la primera ecuacion x, vinculada con el grado de libertad
      %x como dependiente (de ese grado de libertad toma el valor constante), y, segundo, se vincula la
      %ecuacion y con el grado de libertad dependiente con el y.
      %Para interpretar esta restriccion de media de temperatura en todo el dominio igual al valor impuesto en
      %en el archivo de datos, habria que multiplicar ese valor impuesto por el volumen del dominio. Esto
      %puede generar problemas en los problemas multiescala, donde a veces cuando hay poro y no se lo malla,
      %se impone un valor distinto de volumen que no es igual volumen calculado como la suma de los volumenes
      %de los elementos finitos. Por ello se interpreta, que esta restriccion es la interal_dom[Var] = alfa;
      %donde alfa es el valor impuesto en el archivo de datos. Ahora para que esta restriccion tenga sentido
      %de media de Var, alfa debe ser igual a Vol[dom]*ValMedio_Var.
      if ndn==2
      m_RestrCte = [m_ValNodRestr(1,m_IndGdlRestr(1,:));m_ValNodRestr(2,m_IndGdlRestr(2,:))];
      elseif ndn==1
         m_RestrCte = m_ValNodRestr(1,m_IndGdlRestr(1,:));
      else
         error(['Lectura de datos: Condiciones de borde: Determinacion de matrices: Restricciones ',...
            'de media: No esta implementado para que tenga solo dos restricciones de este tipo.'])
      end  
      m_GdlRestrCte = m_GdlDep;
   else
      %En el caso de no tener restricciones se envia matrices vacias, que al concatenar las matrices no
      %influye. Esto se hizo para la homogeneizacion de las deformaciones fuera de la funcion, pero ver si no
      %queda mas ordenado hacerlo dentro de las mismas.
      m_RestrLin = [];
      m_GdlVinc = [];
      m_GdlRestr = [];
      m_RestrCte = [];
      m_GdlRestrCte = [];
   end
   
end