function matlab2gid_mesh(in,xx,e_DatSet,e_VG)
   %
   %  Genracion archivos de post-proces con GID 
   %
   %global filename_res iel
   %iel = e_VG.iel;
   %telem = e_VG.telem;
   %struhyp = e_VG.struhyp;
   %conshyp = e_VG.conshyp;
   % path_file = e_VG.path_file;
   % file = e_VG.file;
   fileCompleto = e_VG.fileCompleto;
   ndime = e_VG.ndime;
   %nElem = e_VG.nElem;
   %nnod = e_VG.nnod;
   ndn = e_VG.ndn;
   %npe = e_VG.npe;
   %eltype = e_VG.eltype;
   %npg = e_VG.npg;
   %xg = e_VG.xg;
   nSet = e_VG.nSet;

   c_elem_type = cell(nSet,1);
   for iSet = 1:nSet

      switch e_DatSet(iSet).e_DatElem.eltype
         case {2,10,32}    % Triangulo de tres nodos 
            c_elem_type{iSet} = 'Triangle';
         case {4,8,20,21,22,23,31,108}  % Cuadrangulos de 4 nodos en desplazamientos
            c_elem_type{iSet} = 'Quadrilateral';
         case 16  % Cuadrangulos de 8 nodos en desplazamientos %AA
            c_elem_type{iSet} = 'Quadrilateral';            %AA
         case 5          % Elemento de barra de 2 nodos en desplazamientos (2D)
            c_elem_type{iSet} = 'Linear';
         case 7          % Elemento hexaedrico de 8 nodos (3D)
            c_elem_type{iSet} = 'Hexahedra';
      %    case 'tet'
      %       elem_type = 'Tetrahedra';
      %       NNode=npe;
         otherwise
            error('PostProceso GiD: Malla: Tipo elemento GiD: Elemento finito no definido.')
      end
      
   end
   
   %Nombre del archivo
   %[~,filename] = fileparts(fileCompleto);
   [dir,filename] = fileparts(fileCompleto);
   %
   if e_VG.isPostResMultFile      
      dir = fullfile(dir,'GiDRes');
      if ~exist(dir,'dir')
         mkdir(dir)
      end
      fileCompleto = fullfile(dir,filename);
   end
   filename_write = [fileCompleto,'.flavia.msh'];
   fid = fopen(filename_write,'wt');
   
   % 1.- Header with 6 free lines
   fprintf(fid,'# ================================================== \n');
   fprintf(fid,'# \n');
   fprintf(fid,'# PostProceso GiD - Archivo de malla \n');
   fprintf(fid,['# Nombre de archivo: ',filename,' \n']);
   fprintf(fid,'# ================================================== \n');
   
   %Para guardar en el GiD la malla, es necesario crear un conjunto MESH por cada tipo de elemento.
   %Como se puede tener distintos SET por ser de materiales distintos, y no por tipo de elemento, no
   %necesit�ndose crear un MESH diferente en el archivo, pero por simplificidad por ahora se crea un
   %MESH por cada SET.
   for iSet = 1:nSet

      npe = e_DatSet(iSet).e_DatElem.npe;
      nElemSet = e_DatSet(iSet).nElem;
      m_Conec = e_DatSet(iSet).conec;
  
      fprintf(fid,['MESH "Set_',num2str(iSet),'" dimension ',num2str(ndime),' Elemtype ',...
         c_elem_type{iSet},' Nnode ',num2str(npe),'\n']);

      % 4.- Coordinates
      fprintf(fid,'Coordinates \n');
      %S�lo es necesario guardar las coordenadas en un solo MESH, por lo que se pone solo en el
      %primero.
      if iSet==1         
         fprintf(fid,'# N�mero_de_nodo Coordenada_x Coordenada_y Coordenada_z \n');
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         fprintf(fid,format,[in,xx(:,1:ndime)]');
      end
      fprintf(fid,'end coordinates \n');
      
      %fprintf(fid,' # we put the same material in the whole MESH, \n');
      fprintf(fid,'Elements \n');

      format = '';
      for iNpe = 1:npe
         format = [format,'Nodo_',num2str(iNpe),' '];  %#ok<AGROW>
      end
      format = ['# Nro_Elemento ',format,' Nro_Set \n'];   %#ok<AGROW>
      fprintf(fid,format);

      format = ['%d',repmat(' %d',1,npe),' %d\n'];
      %La enumeraci�n de los nodos que viene en la conectividad es la local dentro del programa (seg�n el
      %orden en que fue indicado el nodo en la matriz de coordenadas), por lo que se cambia para la numeraci�n
      %de los nodos indicado por el usuario.
      m_Conec(m_Conec~=0) = in(m_Conec(m_Conec~=0));
      fprintf(fid,format,[e_DatSet(iSet).m_NumElem',m_Conec,repmat(iSet,nElemSet,1)]');

      fprintf(fid,'end elements \n');
   
   end

   fclose(fid);

   %%
   % ============================
   %  Archivo file.flavia.res
   % ============================
   %Se guarda el encabezado en el archivo del resultado para el GiD en esta parte ya que
   %matlab2gid_res es llamado en cada paso
   filename_res = [fileCompleto,'.flavia.res'];
   fid_res = fopen(filename_res,'wt');
   fprintf(fid_res,'GiD Post Results File 1.0 \n');
   
   %Como se dijo antes, ser�a m�s �ptimo no escribir lo mismo para cada set, ya que puede varios de
   %ellos que tiene la misma estructura de puntos de gauss, pero por simplificidad por ahora se
   %imprime para cada Set.
   for iSet = 1:nSet
      
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      xg = e_DatElemSet.xg;
      npg = e_DatElemSet.npg;
      
      fprintf(fid,['# Datos de puntos de gauss del Set_',num2str(iSet),'\n']);
      
      %% Definici�n de punto de Gauss
      %Se almacena los puntos de gauss, para no hacerlo en todos los pasos que se plotea los
      %resultados (solo es necesario definirlo una vez).

      %%Punto de Gauss �nico para valor medio de los puntos de gauss
      fprintf(fid_res,['GaussPoints "GP_Unico_Set_',num2str(iSet),'" Elemtype ',...
         c_elem_type{iSet},' "Set_',num2str(iSet),'"\n']);
      fprintf(fid_res,'Number of Gauss Points: 1\n');
      fprintf(fid_res,'Nodes not included\n');
      fprintf(fid_res,'Natural Coordinates: Internal\n');
      fprintf(fid_res,'End GaussPoints\n');

      NPG_ele = npg;
      m_xgElem = xg;
      switch e_DatElemSet.eltype
         case 10    % Tri�ngulo de 3 nodos con discontinuidades fuertes (SDA)
            %Para que el GiD no imprima los resultados de los 2 PG en la misma posici�n, para la
            %impresi�n se desv�a uno del otro levemente, con la misma y 
            %(como se separan se elige arbitrariamente).
            %Ver como hacer correctamente esta impresi�n. Probar que las �ltimas versiones de GiD permite
            %tratar el caso los con PG coincidente.
            %Se cambia la posici�n usando las coordenadas internas del elemento (dimensi�n
            %unitaria).            
            m_xgElem(1,1) = m_xgElem(1,1)*(1-1e-1);
            m_xgElem(2,1) = m_xgElem(2,1)*(1+1e-1);
            %En el caso de imprimir un solo punto de gauss, el regular.
%             m_xgElem = m_xgElem(1,:);
%             NPG_ele = 1;
         case {21,22,23}    % Cuadr�ngulo de 4 nodos con discontinuidades fuertes (SDA)
            NPG_ele = npg-2;
            m_xgElem = m_xgElem(1:NPG_ele,:);
      end
      %Punto de Gauss de los elementos usados (este deber�a ser un loop en todos los Sets)
      fprintf(fid_res,['GaussPoints "GP_Set_',num2str(iSet),'" Elemtype ',...
         c_elem_type{iSet},' "Set_',num2str(iSet),'"\n']);
      fprintf(fid_res,['Number of Gauss Points: ',num2str(NPG_ele),'\n']);
      fprintf(fid_res,'Nodes not included\n');
      %Se utiliza Given en lugar de Internal porque es m�s general para distintos puntos de gauss
      %que se defina, aunque igualmente la cantidad de puntos tiene que ser compatible con lo que
      %permite el GiD para cada elemento y el sistema de coordenadas usada en el mismo.
      fprintf(fid_res,'Natural Coordinates: Given\n');
      fprintf(fid_res,'%f %f\n',m_xgElem');
      fprintf(fid_res,'End GaussPoints\n');
     
      switch e_DatMatSet.conshyp
         %Falta ver como son los fload de los modelos constitutivos 
         case 8
            %% Definici�n de tabla de resultados
            fprintf(fid,['# Tabla de resultados del Set_',num2str(iSet),'\n']);
            %Tabla para los Factores de carga (para verlo en el GID hay que usar Contour Ranges.
            fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de fuerzas centradas"\n']);
            fprintf(fid_res,'- 0: "-1: Descarga el�stica en zona de discontinuidad fuerte"\n');
            fprintf(fid_res,'0 - 1: "0: Zona El�stica"\n');
            fprintf(fid_res,'1 - : "1: Carga en zona de discontinuidad fuerte"\n');
            fprintf(fid_res,'End ResultRangesTable\n');
         case {4,5,10,11,12}
            %% Definici�n de tabla de resultados
            fprintf(fid,['# Tabla de resultados del Set_',num2str(iSet),'\n']);
            %Tabla para los Factores de carga (para verlo en el GID hay que usar Contour Ranges.
            fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de Da�o Isotr�pico"\n']);
            fprintf(fid_res,' - -0.5: "-1: El�stica (con da�o)"\n');
            fprintf(fid_res,'-0.5 - 0.5: "0: El�stica (sin da�o)"\n');
            fprintf(fid_res,'0.5 - 1.5: "1: Da�o"\n');
            fprintf(fid_res,'1.5 - : "2: Sobre l�mite de da�o"\n');
            fprintf(fid_res,'End ResultRangesTable\n');
         case {50,51,53,54,55,60,61}
            f_MallaCUME(e_DatSet(iSet).m_NumElem,e_DatMatSet)
      end
      
   end
   %
   %Cuando se imprime en archivos separados, el Gid tira un peque�o error si un flavia.res tiene los datos del
   %punto gauss y tabla de resultados, sin ning�n resultado, por lo que se imprime un paso 0 sin para evitar
   %este problema.
   if e_VG.isPostResMultFile   
      fprintf(fid_res,'Result "" "Load Analysis" 0 Scalar OnNodes\nValues\nEnd Values\n');
   end
   %%
   fclose(fid_res);
   
   %% Archivo List
   %Inicializaci�n de los archivos lista, que guarda la lista de archivos necesarios para abrir todos los
   %resultados en el caso usarse la opci�n de archivos separados.
   if e_VG.isPostResMultFile
      fileNameLst = [fileCompleto,'.post.lst'];
      fId = fopen(fileNameLst,'wt');
      fprintf(fId,'Merge\n');
      %No es necesario agregar el archivo flavia.res del igual nombre que flavia.msh ya que al leer este
      %�ltimo se lee autom�ticamente el primero.
      fprintf(fId,'%s\n',[filename,'.flavia.msh']);
      fclose(fId);
   end
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if e_VG.protype==1||e_VG.protype==3
       c_elem_type{iSet} = 'Quadrilateral';
   %Nombre del archivo
   %[~,filename] = fileparts(fileCompleto);
   [dir,filename] = fileparts(fileCompleto);
   %
   if e_VG.isPostResMultFile      
      dir = fullfile(dir,'GiDRes');
      if ~exist(dir,'dir')
         mkdir(dir)
      end
      fileCompleto = fullfile(dir,filename);
   end
   filename_writepp = [fileCompleto,'.pp.msh'];
   fidpp = fopen(filename_writepp,'wt');
   
   % 1.- Header with 6 free lines
   fprintf(fidpp,'# ================================================== \n');
   fprintf(fidpp,'# \n');
   fprintf(fidpp,'# PostProceso GiD - Archivo de malla \n');
   fprintf(fidpp,['# Nombre de archivo: ',filename_writepp,' \n']);
   fprintf(fidpp,'# ================================================== \n');
   
   %Para guardar en el GiD la malla, es necesario crear un conjunto MESH por cada tipo de elemento.
   %Como se puede tener distintos SET por ser de materiales distintos, y no por tipo de elemento, no
   %necesitandose crear un MESH diferente en el archivo, pero por simplificidad por ahora se crea un
   %MESH por cada SET.
   for iSet = 1:nSet

      npe = 4;
      nElemSet = e_DatSet(iSet).nElem;
      m_Conec = e_DatSet(iSet).conec(:,1:4);
      omconec = sort(reshape(m_Conec,[],1));
  
      fprintf(fidpp,['MESH "Set_',num2str(iSet),'" dimension ',num2str(ndime),' Elemtype ',...
         c_elem_type{iSet},' Nnode ',num2str(npe),'\n']);

      % 4.- Coordinates
      fprintf(fidpp,'Coordinates \n');
      %Solo es necesario guardar las coordenadas en un solo MESH, por lo que se pone solo en el
      %primero.
      if iSet==1  
          omconec = [];
          for iSet2 = 1:nSet
               m_Conec2 = e_DatSet(iSet2).conec(:,1:4);
               omconec2 = reshape(m_Conec2,[],1);
               omconec = [omconec;omconec2];
          end
          omconec = sort(unique(omconec));
          
         fprintf(fidpp,'# Numero_de_nodo Coordenada_x Coordenada_y Coordenada_z \n');
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         fprintf(fidpp,format,[in(omconec),xx(omconec,1:ndime)]');
      end
      fprintf(fidpp,'end coordinates \n');
      
      %fprintf(fid,' # we put the same material in the whole MESH, \n');
      fprintf(fidpp,'Elements \n');

      format = '';
      for iNpe = 1:npe
         format = [format,'Nodo_',num2str(iNpe),' '];  %#ok<AGROW>
      end
      format = ['# Nro_Elemento ',format,' Nro_Set \n'];   %#ok<AGROW>
      fprintf(fidpp,format);

      format = ['%d',repmat(' %d',1,npe),' %d\n'];
      %La enumeracion de los nodos que viene en la conectividad es la local dentro del programa (segun el
      %orden en que fue indicado el nodo en la matriz de coordenadas), por lo que se cambia para la numeracion
      %de los nodos indicado por el usuario.
      m_Conec(m_Conec~=0) = in(m_Conec(m_Conec~=0));
      fprintf(fidpp,format,[e_DatSet(iSet).m_NumElem',m_Conec,repmat(iSet,nElemSet,1)]');

      fprintf(fidpp,'end elements \n');
   
   end

   fclose(fidpp); 
   
   %%
   % ============================
   %  Archivo file.flavia.res
   % ============================
   %Se guarda el encabezado en el archivo del resultado para el GiD en esta parte ya que
   %matlab2gid_res es llamado en cada paso
   filename_resp = [fileCompleto,'.pp.res'];
   fid_resp = fopen(filename_resp,'wt');
   fprintf(fid_resp,'GiD Post Results File 1.0 \n');
   
   %Como se dijo antes, seria mas optimo no escribir lo mismo para cada set, ya que puede varios de
   %ellos que tiene la misma estructura de puntos de gauss, pero por simplificidad por ahora se
   %imprime para cada Set.
   for iSet = 1:nSet
      
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      xg = e_DatElemSet.xg;
      npg = e_DatElemSet.npg;
      
      fprintf(fid,['# Datos de puntos de gauss del Set_',num2str(iSet),'\n']);
      
      %% Definicion de punto de Gauss
      %Se almacena los puntos de gauss, para no hacerlo en todos los pasos que se plotea los
      %resultados (solo es necesario definirlo una vez).

      %%Punto de Gauss unico para valor medio de los puntos de gauss
      fprintf(fid_resp,['GaussPoints "GP_Unico_Set_',num2str(iSet),'" Elemtype ',...
         c_elem_type{iSet},' "Set_',num2str(iSet),'"\n']);
      fprintf(fid_resp,'Number of Gauss Points: 1\n');
      fprintf(fid_resp,'Nodes not included\n');
      fprintf(fid_resp,'Natural Coordinates: Internal\n');
      fprintf(fid_resp,'End GaussPoints\n');

      NPG_ele = npg;
      m_xgElem = xg;
      %Punto de Gauss de los elementos usados (este deberia ser un loop en todos los Sets)
      fprintf(fid_resp,['GaussPoints "GP_Set_',num2str(iSet),'" Elemtype ',...
         c_elem_type{iSet},' "Set_',num2str(iSet),'"\n']);
      fprintf(fid_resp,['Number of Gauss Points: ',num2str(NPG_ele),'\n']);
      fprintf(fid_resp,'Nodes not included\n');
      %Se utiliza Given en lugar de Internal porque es mas general para distintos puntos de gauss
      %que se defina, aunque igualmente la cantidad de puntos tiene que ser compatible con lo que
      %permite el GiD para cada elemento y el sistema de coordenadas usada en el mismo.
      fprintf(fid_resp,'Natural Coordinates: Given\n');
      fprintf(fid_resp,'%f %f\n',m_xgElem');
      fprintf(fid_resp,'End GaussPoints\n');
     
      switch e_DatMatSet.conshyp
          case {50,51,53,54,55,60,61}
            f_MallaCUME(e_DatSet(iSet).m_NumElem,e_DatMatSet)
      end
      
   end
   %
   %Cuando se imprime en archivos separados, el Gid tira un pequeño error si un flavia.res tiene los datos del
   %punto gauss y tabla de resultados, sin ningun resultado, por lo que se imprime un paso 0 sin para evitar
   %este problema.
   if e_VG.isPostResMultFile   
      fprintf(fid_resp,'Result "" "Load Analysis" 0 Scalar OnNodes\nValues\nEnd Values\n');
   end
   %%
   fclose(fid_resp);
   end
%    %% Archivo List
%    %Inicializacion de los archivos lista, que guarda la lista de archivos necesarios para abrir todos los
%    %resultados en el caso usarse la opcion de archivos separados.
%    if e_VG.isPostResMultFile
%       fileNameLst = [fileCompleto,'.post.lst'];
%       fId = fopen(fileNameLst,'wt');
%       fprintf(fId,'Merge\n');
%       %No es necesario agregar el archivo flavia.res del igual nombre que flavia.msh ya que al leer este
%       %ultimo se lee automaticamente el primero.
%       fprintf(fId,'%s\n',[filename,'.flavia.msh']);
%       fclose(fId);
%    end 
end