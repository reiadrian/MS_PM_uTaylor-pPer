function f_InicArchDat(in,m_SetElem,e_DatSet,e_VG)

   nSet = e_VG.nSet;
   fileCompleto = e_VG.fileCompleto;
   
   %% Archivo que indica que paso del esquema del newton que no convergi�.
   nomArch = [fileCompleto,'.pasNoConv'];
   %Se borra para que cuando no converja se cree el archivo y indique que no convergi� alg�n paso.
   %Por ahora no se crea estos archivos para el problema multiescala, ya que si son mucho los elementos de
   %este tipo genera muchos archivos que hacen lento el directorio y tarda mucho. Tambi�n va ser un problema
   %en el caso de que se use en un cluster y no encuentre el archivo donde escribir.
%    if ~e_VG.esME
%       fId = fopen(nomArch,'wt');
%       fprintf(fId,'#Pasos en que el esquema de newton no convergi�.\n');
%       fclose(fId);
%    end
   %Para que el caso multiescala se imprima un solo archivo de pasos de tiempo para todos los PG macro, se
   %se usa un nombre �nico.
   if e_VG.esME
      nomArch = [e_VG.fileCompletoOrig,'.pasNoConv'];
   end
   fId = fopen(nomArch,'wt');
   fprintf(fId,'#Pasos en que el esquema de newton no convergi�.\n');
   fclose(fId);
   
   %% Archivos seg�n el elemento y el modelo constitutivo
   for iSet = 1:nSet
      
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      m_NumElem = e_DatSet(iSet).m_NumElem;      
      nElem = e_DatSet(iSet).nElem;
      npg = e_DatElemSet.npg;
            
      for iElem = 1:nElem
         
         eltype = e_DatElemSet.eltype;
         conshyp = e_DatMatSet.conshyp;
      
         for iPG = 1:npg
            
            %Seg�n el tipo de elemento se inicializa los archivos de an�lisis de bifurcaci�n.
            switch eltype
               case {2,4,10,22,23,108}  % Tri�ngulo de 3 nodos con discontinuidades fuertes (SDA)
                  if iPG==1
                     %An�lisis de bifurcaci�n
                     %Se asume que este elemento tiene dos puntos de gauss, en la misma
                     %posici�n, por lo que se imprime el an�lisis de bifurcaci�n del punto regular.
%                      nombrArchBif = [fileCompleto,'_E',int2str(m_NumElem(iElem)),'.analisisBif']; 
%                      fId = fopen(nombrArchBif,'wt');
%                      fprintf(fId,'#An�lisis de bifurcaci�n.\n');
%                      fclose(fId);
                  end
                  %Para tener un solo archivo para el todo set.
                  if iElem==1&&iPG==1
                     %An�lisis de bifurcaci�n
                     %Este se asume que este elemento tiene dos puntos de gauss, en la misma
                     %posici�n, por lo que se imprime el an�lisis de bifurcaci�n del punto regular.
                     %Notar que si se tiene varios sets con el mismo modelo constitutivo se utiliza
                     %el mismo archivo (ver si no cambiar el nombre por set)
                     nomArch = [fileCompleto,'.elemBif'];
                     fId = fopen(nomArch,'wt');
                     fprintf(fId,'#Elementos que bifurcaron y paso que lo hicieron.\n');
                     fprintf(fId,['#| Paso (tiempo) | Elemento | �ngulos de la normal ',...
                        'usado (�ngulos del an�lisis de bifurcaci�n)\n']);
                     fclose(fId);
                  end
            end
            
            %Seg�n el tipo modelo constitutivo se inicializa los archivos de datos
            switch conshyp
               case {50,51,60,61.62} %Modelo multiescala cl�sico y cohesivo
                  e_VGMicro = e_DatMatSet.e_VG;
                  e_VGMicro.fileCompletoOrig = e_VGMicro.fileCompleto;
                  e_VGMicro.fileCompleto = [e_VGMicro.fileCompleto,'_EM',...
                     int2str(m_NumElem(iElem)),'PGM',int2str(iPG)];
                  e_VGMicro.esME = true;
                  f_InicArchDat(e_DatMatSet.in,e_DatMatSet.m_SetElem,e_DatMatSet.e_DatSet,e_VGMicro)
                  %Se imprime el an�lisis de bifurcaci�n para el punto de gauss regular.
                  if (eltype==10 ||eltype==22||eltype==23) &&iPG==1&&any(e_DatMatSet.m_ElemPGImpr(...
                        1,e_DatMatSet.m_ElemPGImpr(2,:)==iPG)==iElem)
                     %An�lisis de bifurcaci�n
                     nombrArchBif = [fileCompleto,'_E',int2str(m_NumElem(iElem)),'.analisisBif']; 
                     fId = fopen(nombrArchBif,'wt');
                     fprintf(fId,'#An�lisis de bifurcaci�n.\n');
                     fclose(fId);
                  end
            end
            
         end
      
      end
      
   end
   
   %% Inicializaci�n del del archivo de datos y impresi�n del encabezado
%    c_NomDat = {'T','Dx','Dy','Dz','Fx','Fy','Fz','Bx','By','Bz','Tx','Ty','Tz',...
%       'Exx','Eyy','Ezz','Exy','Eyx','Txx','Tyy','Tzz','Txy','Tyx','Efxx','Efyy','Efzz','Efxy',...
%       'Da'};
%    c_NomDat = {'T','Dx','Dy','Dz','Fx','Fy','Fz','Bx','By','Bz','Tx','Ty','Tz',...
%       'Exx','Eyy','Ezz','Exy','Eyx','Txx','Tyy','Tzz','Txy','Tyx','Efxx','Efyy','Efzz','Efxy',...
%       'Da','DisGLO','Ehxx','Ehyy','Ehzz','Ehxy','Shxx','Shyy','Shzz','Shxy','Snn','Snt','p','q'}; %JLM
   c_NomDat = {'T','Dx','Dy','Dz','Fx','Fy','Fz','Bx','By','Bz','Tx','Ty','Tz',...
      'Exx','Eyy','Ezz','Exy','Eyx','Txx','Tyy','Tzz','Txy','Tyx','TExx','TEyy','TEzz','TExy','TEyx',... %AA: add TE
      'TTxx','TTyy','TTzz','TTxy','TTyx','Porepress','Efxx','Efyy','Efzz','Efxy',... %AA: add TT and Porepress
      'Da','DisGLO','Ehxx','Ehyy','Ehzz','Ehxy','Shxx','Shyy','Shzz','Shxy','Snn','Snt','p','q','logT',...
      'Vx','Vy','Vx_n+theta','Vy_n+theta','chi'}; %JLM %AA:logT

  c_TextTipoEjeGraf = {'X','Y'};

   %En m_DatGrafXY se tiene organizados los datos de la siguiente manera:
   %{'Nx','Ny','Ex','Ey','PGx','PGy','X','Y'}
   m_DatGrafXY = e_VG.m_DatGrafXY;      
      
   if ~isempty(m_DatGrafXY)
      nGraf = size(m_DatGrafXY,1);
      c_TextoEjes = cell(2,1);
      for iGraf = 1:nGraf         
         for iEje = 1:2
            tipoDat = c_NomDat{m_DatGrafXY(iGraf,iEje+6)};
            switch tipoDat
               case 'T'                     
                  c_TextoEjes{iEje} = 'Tiempo';
               case 'Dx'
                  nodo = in(m_DatGrafXY(iGraf,iEje));
                  f_VerifInfoGraf(nodo,'nodo',c_TextTipoEjeGraf{iEje},iGraf);
                  c_TextoEjes{iEje} = ['Componente_X_del_desplazamiento_del_nodo ',int2str(nodo)];
               case 'Dy'
                  nodo = in(m_DatGrafXY(iGraf,iEje));
                  f_VerifInfoGraf(nodo,'nodo',c_TextTipoEjeGraf{iEje},iGraf);
                  c_TextoEjes{iEje} = ['Componente_Y_del_desplazamiento_del_nodo ',int2str(nodo)];
               case 'Fx'
                  nodo = in(m_DatGrafXY(iGraf,iEje));
                  f_VerifInfoGraf(nodo,'nodo',c_TextTipoEjeGraf{iEje},iGraf);
                  c_TextoEjes{iEje} = ['Componente_X_de_la_fuerza_interna_del_nodo ',int2str(nodo)];
               case 'Fy'
                  nodo = in(m_DatGrafXY(iGraf,iEje));
                  f_VerifInfoGraf(nodo,'nodo',c_TextTipoEjeGraf{iEje},iGraf);
                  c_TextoEjes{iEje} = ['Componente_Y_de_la_fuerza_interna_del_nodo ',int2str(nodo)];
               case 'Bx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Componente_X_del_salto_del_elemento ',int2str(elem)];
               case 'By'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Componente_Y_del_salto_del_elemento ',int2str(elem)];
               case 'Tx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Componente_X_de_la_tracci�n_del_elemento ',int2str(elem)];
               case 'Ty'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Componente_Y_de_la_tracci�n_del_elemento ',int2str(elem)];
               case 'Exx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_Exx_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Eyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_Eyy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Ezz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_Ezz_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Exy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_Exy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss',int2str(pg)];
               case 'Eyx'     %Para el caso de tensores no sim�tricos (LD por ejemplo)
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_Eyx_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];       
               case 'Efxx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_fluctuante_Exx_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Efyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss ',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_fluctuante_Eyy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Efzz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_fluctuante_Ezz_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Efxy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);                  
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Deformacion_fluctuante_Exy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Txx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_Txx_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Tyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_Tyy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Tzz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_Tzz_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Txy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_Txy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Tyx'     %Para el caso de tensores no sim�tricos (LD por ejemplo)
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_Tyx_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
                % AA Tensiones efectivas
                case 'TExx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_efectiva_TExx_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'TEyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_efectiva_TEyy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'TEzz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_efectiva_TEzz_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'TExy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_efectiva_TExy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
                % AA
                % AA Tensiones totales
                case 'TTExx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_total_TTxx_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'TTyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_total_TTyy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'TTzz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_total_TTzz_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'TTxy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Tension_total_TTxy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
                % AA
                % AA Poropresiones
                case 'Porepress'
                  nodo = in(m_DatGrafXY(iGraf,iEje));
                  f_VerifInfoGraf(nodo,'nodo',c_TextTipoEjeGraf{iEje},iGraf);
                  c_TextoEjes{iEje} = ['Propresion_del_nodo ',int2str(nodo)];
               % AA
               % AA Logaritmo del tiempo
               case 'logT'                     
                  c_TextoEjes{iEje} = 'logT';
               % AA Tensiones efectivas
                case 'Vx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['VelocidadFiltracion_Vx_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Vy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['VelocidadFiltracion_Vy_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
             case 'Vx_n+theta'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['VelocidadFiltracion_Vx_n+theta_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Vy_n+theta'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['VelocidadFiltracion_Vy_n+theta_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'chi'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Masa_fluido_chi_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'Da'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  f_VerifInfoGraf(elem,'elemento',c_TextTipoEjeGraf{iEje},iGraf);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  f_VerifInfoGraf(pg,'punto de gauss',c_TextTipoEjeGraf{iEje},iGraf);
                  set = m_SetElem(elem);
                  elem = e_DatSet(set).m_NumElem(e_DatSet(set).m_IndElemSet==elem);
                  c_TextoEjes{iEje} = ['Variable_de_da�o_del_elemento ',int2str(elem),...
                     ' y_del_punto_de_gauss ',int2str(pg)];
               case 'DisGLO'
                  c_TextoEjes{iEje} = ('Disipacion_estructural_del_test');
                case 'Ehxx'
                  c_TextoEjes{iEje} = ('Homogenized_Strain_(component_xx)');
                case 'Ehyy'
                    c_TextoEjes{iEje} = ('Homogenized_Strain_(component_yy)');
                case 'Ehzz'
                    c_TextoEjes{iEje} = ('Homogenized_Strain_(component_zz)');
                case 'Ehxy'
                    c_TextoEjes{iEje} = ('Homogenized_Strain_(component_xy)');
                case 'Shxx'
                    c_TextoEjes{iEje} = ('Homogenized_Stress_(component_xx)');
                case 'Shyy'
                    c_TextoEjes{iEje} = ('Homogenized_Stress_(component_yy)');
                case 'Shzz'
                    c_TextoEjes{iEje} = ('Homogenized_Stress_(component_zz)');
                case 'Shxy'
                    c_TextoEjes{iEje} = ('Homogenized_Stress_(component_xy)');
                case 'Snn'
                    c_TextoEjes{iEje} = ('Homogenized_Stress_(component_nn)');
                case 'Snt'
                    c_TextoEjes{iEje} = ('Homogenized_Stress_(component_nt)');
                case 'p'
                    c_TextoEjes{iEje} = ('Mean_spherical_stress');
                case 'q'
                    c_TextoEjes{iEje} = ('J2');
               otherwise
                  error('Archivos de datos: Inicializaci�n: No est� definido este tipo de dato.')
            end
         end
%          fId = fopen([fileCompleto,'.cur',num2str(iGraf,'%03d')],'wt'); %AA
         fId = fopen([fileCompleto,num2str(iGraf,'%03d'),'.dat'],'wt'); %AA
         %Se usa el s�mbolo de comentario est�ndar del GNUPlot
%          fprintf(fId,['###',c_TextoEjes{1},' Vs ',c_TextoEjes{2},'###''\n']); %AA
         if iGraf==1
         fprintf(fId,['Numero_de_graficos: ',num2str(nGraf,'%d'),'\n']); %AA
         end
         if e_VG.protype==0
         fprintf(fId,['Steps ',num2str(e_VG.np,'%d'),'\n']); %AA
         elseif e_VG.protype==1
         fprintf(fId,['Steps ',num2str(e_VG.np+1,'%d'),'\n']); %AA  
         end
         fprintf(fId,[c_TextoEjes{1},'\n']); %AA
         fprintf(fId,['VS','\n']); %AA
         fprintf(fId,[c_TextoEjes{2},'\n']); %AA
         fclose(fId);            
      end         
   end

end

function f_VerifInfoGraf(ubic,textTipoUbic,textTipoEje,iGraf)
   
   if isnan(ubic)
      error(['Archivos de datos: Inicializaci�n: Se debe ingresar el n�mero ',...
         'de %s para los datos del eje %s de la gr�fica %d.'],textTipoUbic,textTipoEje,iGraf)
   end
   
end