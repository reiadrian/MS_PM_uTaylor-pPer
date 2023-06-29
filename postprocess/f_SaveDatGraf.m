%function f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst,e_VarAux,e_DatSet,m_SetElem,e_VG) %AA
function f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst,e_VarAux,e_DatSet,m_SetElem,...
             sigmaHomog,epsilon_Macro,e_VG)
   % Impresi�n de los resultados de los datos pedidos

   %Celda con los nombres de los datos
   
%    c_NomDat = {'T','Dx','Dy','Dz','Fx','Fy','Fz','Bx','By','Bz','Tx','Ty','Tz',...
%       'Exx','Eyy','Ezz','Exy','Eyx','Txx','Tyy','Tzz','Txy','Tyx','Efxx','Efyy','Efzz','Efxy',...
%       'Da'};
%   c_NomDat = {'T','Dx','Dy','Dz','Fx','Fy','Fz','Bx','By','Bz','Tx','Ty','Tz',...
%       'Exx','Eyy','Ezz','Exy','Eyx','Txx','Tyy','Tzz','Txy','Tyx','TExx','TEyy','TEzz','TExy','TEyx',... %AA: add TE
%       'TTxx','TTyy','TTzz','TTxy','TTyx','Efxx','Efyy','Efzz','Efxy',... %AA: add TT
%       'Da','DisGLO','Ehxx','Ehyy','Ehzz','Ehxy','Shxx','Shyy','Shzz','Shxy','Snn','Snt','p','q'}; %JLM

    c_NomDat = {'T','Dx','Dy','Dz','Fx','Fy','Fz','Bx','By','Bz','Tx','Ty','Tz',...
      'Exx','Eyy','Ezz','Exy','Eyx','Txx','Tyy','Tzz','Txy','Tyx','TExx','TEyy','TEzz','TExy','TEyx',... %AA: add TE
      'TTxx','TTyy','TTzz','TTxy','TTyx','Porepress','Efxx','Efyy','Efzz','Efxy',... %AA: add TT and Porepress
      'Da','DisGLO','Ehxx','Ehyy','Ehzz','Ehxy','Shxx','Shyy','Shzz','Shxy','Snn','Snt','p','q','logT',...
      'Vx','Vy','Vx_n+theta','Vy_n+theta','chi'}; %JLM %AA:logT

   %En m_DatGrafXY se tiene organizados los datos de la siguiente manera:
   %{'Nx','Ny','Ex','Ey','PGx','PGy','X','Y'}
   m_DatGrafXY = e_VG.m_DatGrafXY;
   fileCompleto = e_VG.fileCompleto;
      
   if ~isempty(m_DatGrafXY)
      nGraf = size(m_DatGrafXY,1);
      m_Dat = zeros(1,2);
      for iGraf = 1:nGraf
         for iEje = 1:2
            tipoDat = c_NomDat{m_DatGrafXY(iGraf,iEje+6)};            
            switch tipoDat
               case 'T'                     
                  m_Dat(iEje) = e_VG.Dtime*e_VG.istep;
               case 'Dx'
                  iGdl = 1;
                  nodo = m_DatGrafXY(iGraf,iEje);
                  gdl = e_VG.ndn*nodo-(e_VG.ndn-iGdl);
                  m_Dat(iEje) = u(gdl);
               case 'Dy'
                  iGdl = 2;
                  nodo = m_DatGrafXY(iGraf,iEje);
                  gdl = e_VG.ndn*nodo-(e_VG.ndn-iGdl);
                  m_Dat(iEje) = u(gdl);
               case 'Fx'
                  iGdl = 1;
                  %Tener cuidado que si hay variables condensadas (como el salto), la Fint del nodo
                  %no corresponde al valor de la reacci�n, sino que a una Fint*, es decir la Fint
                  %correspondiente de la condesaci�n.
                  nodo = m_DatGrafXY(iGraf,iEje);
                  gdl = e_VG.ndn*nodo-(e_VG.ndn-iGdl);
                  m_Dat(iEje) = Fint(gdl);
               case 'Fy'
                  iGdl = 2;
                  %Tener cuidado que si hay variables condensadas (como el salto), la Fint del nodo
                  %no corresponde al valor de la reacci�n (nodo con restricci�n), sino que a una 
                  %Fint*, es decir la Fint correspondiente a la condesaci�n.
                  nodo = m_DatGrafXY(iGraf,iEje);
                  gdl = e_VG.ndn*nodo-(e_VG.ndn-iGdl);
                  m_Dat(iEje) = Fint(gdl);
               case 'Bx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  set = m_SetElem(elem);
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eltype = e_DatSet(set).e_DatElem.eltype;
                  switch eltype
                     case 10
                        m_Beta = c_GdlCond{set,1};
                        m_Dat(iEje) = m_Beta(1,m_IndElemSet);
                     case 23
                        m_Beta = c_GdlCond{set,1};
                        m_Dat(iEje) = m_Beta(1,m_IndElemSet);
                     otherwise
                        error(['Archivos de datos: Impresi�n: El salto Bx no est� definido',...
                           'para elementos tipo %d.'],eltype)
                  end
               case 'By'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  set = m_SetElem(elem);
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eltype = e_DatSet(set).e_DatElem.eltype;
                  switch eltype
                     case 10
                        m_Beta = c_GdlCond{set,1};
                        m_Dat(iEje) = m_Beta(2,m_IndElemSet);
                     case 23
                        m_Beta = c_GdlCond{set,1};
                        m_Dat(iEje) = m_Beta(2,m_IndElemSet);
                     otherwise
                        error(['Archivos de datos: Impresi�n: El salto By no est� definido',...
                           'para elementos tipo %d.'],eltype)
                  end
               case 'Tx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  set = m_SetElem(elem);
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eltype = e_DatSet(set).e_DatElem.eltype;
                  ndime = e_VG.ndime;
                  switch eltype
                     case 10
                        %m_Dat(iEje) = e_VarAux(set).VarAuxElem(19,m_IndElemSet);
                        m_Dat(iEje) = e_VarEst(set).VarHistElem(ndime+1,m_IndElemSet);
                     case 23
                        p_Traccion = e_DatSet(set).e_DatElem.pointersVAE.p_Traccion;
                        m_Dat(iEje) = e_VarAux(set).VarAuxElem(p_Traccion(1),m_IndElemSet);
                     otherwise
                        error(['Archivos de datos: Impresi�n: La tracci�n Tx no est� definida',...
                           'para elementos tipo %d.'],eltype)
                  end
               case 'Ty'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  set = m_SetElem(elem);
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eltype = e_DatSet(set).e_DatElem.eltype;
                  ndime = e_VG.ndime;
                  switch eltype
                     case 10
                        %m_Dat(iEje) = e_VarAux(set).VarAuxElem(20,m_IndElemSet);
                        m_Dat(iEje) = e_VarEst(set).VarHistElem(2*ndime,m_IndElemSet);
                     case 23
                        p_Traccion = e_DatSet(set).e_DatElem.pointersVAE.p_Traccion;
                        m_Dat(iEje) = e_VarAux(set).VarAuxElem(p_Traccion(2),m_IndElemSet);
                     otherwise
                        error(['Archivos de datos: Impresi�n: La tracci�n Ty no est� definida',...
                           'para elementos tipo %d.'],eltype)
                  end
               case 'Exx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps = reshape(e_VarEst(set).eps,ntens,npg,[]);
                  m_Dat(iEje) = eps(1,pg,m_IndElemSet);
               case 'Eyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps = reshape(e_VarEst(set).eps,ntens,npg,[]);
                  m_Dat(iEje) = eps(2,pg,m_IndElemSet);
               case 'Ezz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps = reshape(e_VarEst(set).eps,ntens,npg,[]);
                  m_Dat(iEje) = eps(3,pg,m_IndElemSet);
               case 'Exy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps = reshape(e_VarEst(set).eps,ntens,npg,[]);
                  m_Dat(iEje) = eps(4,pg,m_IndElemSet);
               case 'Eyx'   %Para el caso de tensores no sim�tricos (LD por ejemplo)
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps = reshape(e_VarEst(set).eps,ntens,npg,[]);
                  m_Dat(iEje) = eps(5,pg,m_IndElemSet); 
               case 'Efxx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps_fluct = reshape(e_VarEst(set).eps_fluct,ntens,npg,[]);
                  m_Dat(iEje) = eps_fluct(1,pg,m_IndElemSet);
               case 'Efyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps_fluct = reshape(e_VarEst(set).eps_fluct,ntens,npg,[]);
                  m_Dat(iEje) = eps_fluct(2,pg,m_IndElemSet);
               case 'Efzz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps_fluct = reshape(e_VarEst(set).eps_fluct,ntens,npg,[]);
                  m_Dat(iEje) = eps_fluct(3,pg,m_IndElemSet);
               case 'Efxy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  eps_fluct = reshape(e_VarEst(set).eps_fluct,ntens,npg,[]);
                  m_Dat(iEje) = eps_fluct(4,pg,m_IndElemSet);
               case 'Txx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigma,ntens,npg,[]);
                  m_Dat(iEje) = tens(1,pg,m_IndElemSet);
               case 'Tyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigma,ntens,npg,[]);
                  m_Dat(iEje) = tens(2,pg,m_IndElemSet);
               case 'Tzz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigma,ntens,npg,[]);
                  m_Dat(iEje) = tens(3,pg,m_IndElemSet);
               case 'Txy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigma,ntens,npg,[]);
                  m_Dat(iEje) = tens(4,pg,m_IndElemSet);
               case 'Tyx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigma,ntens,npg,[]);
                  m_Dat(iEje) = tens(5,pg,m_IndElemSet);
               %AA Tensiones efectivas
               case 'TExx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigmaE,ntens,npg,[]);
                  m_Dat(iEje) = tens(1,pg,m_IndElemSet);
               case 'TEyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigmaE,ntens,npg,[]);
                  m_Dat(iEje) = tens(2,pg,m_IndElemSet);
               case 'TEzz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigmaE,ntens,npg,[]);
                  m_Dat(iEje) = tens(3,pg,m_IndElemSet);
               case 'TExy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigmaE,ntens,npg,[]);
                  m_Dat(iEje) = tens(4,pg,m_IndElemSet);
%                case 'TEyx' %NO EXISTE. ES SIM�TRICO
%                   elem = m_DatGrafXY(iGraf,iEje+2);
%                   pg = m_DatGrafXY(iGraf,iEje+4);
%                   set = m_SetElem(elem);
%                   ntens = e_VG.ntens;
%                   npg = e_DatSet(set).e_DatElem.npg;
%                   m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
%                   tens = reshape(e_VarEst(set).sigmaE,ntens,npg,[]);
%                   m_Dat(iEje) = tens(5,pg,m_IndElemSet);
               %AA   
               %AA Tensiones totales
               case 'TTxx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigmaT,ntens,npg,[]);
                  m_Dat(iEje) = tens(1,pg,m_IndElemSet);
               case 'TTyy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigmaT,ntens,npg,[]);
                  m_Dat(iEje) = tens(2,pg,m_IndElemSet);
               case 'TTzz'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigmaT,ntens,npg,[]);
                  m_Dat(iEje) = tens(3,pg,m_IndElemSet);
               case 'TTxy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigmaT,ntens,npg,[]);
                  m_Dat(iEje) = tens(4,pg,m_IndElemSet);
%                case 'TTyx' %NO EXISTE. ES SIM�TRICO
%                   elem = m_DatGrafXY(iGraf,iEje+2);
%                   pg = m_DatGrafXY(iGraf,iEje+4);
%                   set = m_SetElem(elem);
%                   ntens = e_VG.ntens;
%                   npg = e_DatSet(set).e_DatElem.npg;
%                   m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
%                   tens = reshape(e_VarEst(set).sigmaT,ntens,npg,[]);
%                   m_Dat(iEje) = tens(5,pg,m_IndElemSet);
               %AA Poropresiones
               case 'Porepress'
                  iGdl = 3;
                  nodo = m_DatGrafXY(iGraf,iEje);
                  gdl = e_VG.ndn*nodo-(e_VG.ndn-iGdl);
                  m_Dat(iEje) = u(gdl); 
               %AA
               %AA Logaritmo del tiempo
               case 'logT'
                   logT = f_logT(e_VG,iGraf);
                   m_Dat(iEje) = logT;
               %AA
               %AA Velocidades de filtracion en el tiempo n+1
               case 'Vx'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  nvel = 2;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  vel = reshape(e_VarEst(set).velflu_total,nvel,npg,[]);
                  m_Dat(iEje) = vel(1,pg,m_IndElemSet);
               case 'Vy'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  nvel = 2;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  vel = reshape(e_VarEst(set).velflu_total,nvel,npg,[]);
                  m_Dat(iEje) = vel(2,pg,m_IndElemSet);
               %AA Velocidades de filtracion en el tiempo n+1
               case 'Vx_n+theta'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  nvel = 2;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  vel = reshape(e_VarEst(set).velflu,nvel,npg,[]);
                  m_Dat(iEje) = vel(1,pg,m_IndElemSet);
               case 'Vy_n+theta'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  nvel = 2;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  vel = reshape(e_VarEst(set).velflu,nvel,npg,[]);
                  m_Dat(iEje) = vel(2,pg,m_IndElemSet);
               case 'chi'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  nchi = 1;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  chiX = reshape(e_VarEst(set).mflu,nchi,npg,[]);
                  m_Dat(iEje) = chiX(1,pg,m_IndElemSet);
               %AA
               case 'p'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigma(:,elem),ntens,npg,[]);
                  m_Dat(iEje) = (tens(1,1)+ tens(2,1)+tens(3,1))/3; %tens(5,pg,m_IndElemSet);
               case 'q'
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  ntens = e_VG.ntens;
                  npg = e_DatSet(set).e_DatElem.npg;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  tens = reshape(e_VarEst(set).sigma(:,elem),ntens,npg,[]);
                  J2 = 1/6*((tens(1,1)-tens(2,1))^2 + (tens(2,1)-tens(3,1))^2 + (tens(3,1)-tens(1,1))^2 + ...
                      6*(tens(4,1)^2)) ; 
                  m_Dat(iEje) = J2 ;

%                case {'TxIDL','TyIDL','vTIDL'}
%                   %Tracciones instant�neas en el dominio localizado micro.
%                   elem = m_DatGrafXY(iGraf,iEje+2);
%                   pg = m_DatGrafXY(iGraf,iEje+4);
%                   set = m_SetElem(elem);
%                   %Datos materiales micro
%                   e_DatMat = e_DatSet(set).e_DatMat;
%                   e_DatSetMicro = e_DatMat.e_DatSet;
%                   %Se asume que el modelo cohesivo se tiene un dato por PG en las variables auxiliares y
%                   %hist�ricas, como se maneja estructuras.
%                   %Se utiliza el dominio localizado de c�lculo que es el var�a permanentemente.
%                   m_ElemLocCalc = e_VarAux(set).VarAuxGP(1,pg,elem).m_ElemLocCalc;
%                   omegaMicroL = sum(e_DatMat.m_VolElem(m_ElemLocCalc));
%                   %Se obtiene la normal macro, se asume que es el elemento T1-SDA.
%                   m_Normal = reshape(e_VarAux(set).VarAuxElem(2:9,elem),ntens,ndime); 
%                   %Se utiliza las tensiones de los PG, notar que esto no es correcto para el caso del elemento
%                   %Q1 mixto con inyecci�n de deformaci�n (en este caso hay que usar las tensiones
%                   %estabilizadas), a nivel micro.
%                   c_DetJTLoc = arrayfun(@(x)bsxfun(@times,x.m_DetJT,m_ElemLocCalc(x.m_IndElemSet)),...
%                      e_DatSet(set),'UniformOutput',false);
%                   m_SigmaHomog = f_HomogArea({e_VarEst(set).sigma},ntens,omegaMicroL,c_DetJTLoc,...
%                      e_DatSetMicro,e_VG);
%                   switch tipoDat
%                      case 'TxIDL'
%                         m_Dat(iEje) = m_Normal(:,1)'*m_SigmaHomog;
%                      case 'TyIDL'
%                         m_Dat(iEje) = m_Normal(:,2)'*m_SigmaHomog;
%                      case 'vTIDL'
%                         m_Dat(iEje) = m_Normal'*m_SigmaHomog;
%                   end
%                case 'TyIDL'
%                   %Tracciones instant�neas en el dominio localizado micro, componente y.
%                   elem = m_DatGrafXY(iGraf,iEje+2);
%                   pg = m_DatGrafXY(iGraf,iEje+4);
%                   set = m_SetElem(elem);
%                   %Datos materiales micro
%                   e_DatMat = e_DatSet(set).e_DatMat;
%                   e_DatSetMicro = e_DatMat.e_DatSet;
%                   %Se asume que el modelo cohesivo se tiene un dato por PG en las variables auxiliares y
%                   %hist�ricas, como se maneja estructuras.
%                   %Se utiliza el dominio localizado de c�lculo que es el var�a permanentemente.
%                   m_ElemLocCalc = e_VarAux(set).VarAuxGP(1,pg,elem).m_ElemLocCalc;
%                   omegaMicroL = sum(e_DatMat.m_VolElem(m_ElemLocCalc));
%                   %Se obtiene la normal macro, se asume que es el elemento T1-SDA.
%                   m_Normal = reshape(e_VarAux(set).VarAuxElem(2:9,elem),ntens,ndime); 
%                   %Se utiliza las tensiones de los PG, notar que esto no es correcto para el caso del elemento
%                   %Q1 mixto con inyecci�n de deformaci�n (en este caso hay que usar las tensiones
%                   %estabilizadas) a nivel micro.
%                   c_DetJTLoc = arrayfun(@(x)bsxfun(@times,x.m_DetJT,m_ElemLocCalc(x.m_IndElemSet)),...
%                      e_DatSet(set),'UniformOutput',false);
%                   m_SigmaHomog = f_HomogArea({e_VarEst(set).sigma},ntens,omegaMicroL,c_DetJTLoc,...
%                      e_DatSetMicro,e_VG);
%                   m_Traccy = m_Normal(:,2)'*m_SigmaHomog;
%                   m_Dat(iEje) = m_Traccy;
               case 'Da'
                  %Para los modelos constitutivos con da�o se imprime la variable de da�o.
                  elem = m_DatGrafXY(iGraf,iEje+2);
                  pg = m_DatGrafXY(iGraf,iEje+4);
                  set = m_SetElem(elem);
                  npg = e_DatSet(set).e_DatElem.npg;
                  conshyp = e_DatSet(set).e_DatMat.conshyp;
                  m_IndElemSet = e_DatSet(set).m_IndElemSet==elem;
                  siavarpg = e_DatSet(set).e_DatMat.siavarpg;
                  switch conshyp
                     case {10,11,12,13}                        
                        m_dam = reshape(e_VarAux(set).VarAuxGP,siavarpg,npg,[]);
                        m_Dat(iEje) = m_dam(1,pg,m_IndElemSet);
                     otherwise
                        error(['Archivos de datos: Impresi�n: La variable de da�o no est� definida ',...
                           'para el modelos constitutivo tipo %d.'],conshyp)
                  end

               
                case 'DisGLO'
                    % Disipacion estructural del test
                    DisGlob = 0;
                    for iSet=1:e_VG.nSet
                        if e_DatSet(iSet).e_DatElem.eltype == 21
                            i_VarDissipation  = e_DatSet(iSet).e_DatElem.pointersVHE.i_VarDissipation ;
                            m_VarHistElemNew  = e_VarEst(iSet).VarHistElem;
                            m_dissipation_new = m_VarHistElemNew (i_VarDissipation,:)  ;
                            DisGlob = DisGlob + sum(m_dissipation_new) ;
                        end
                    end
                    m_Dat(iEje) = DisGlob;
                case 'Ehxx'
                    m_Dat(iEje) = epsilon_Macro(1);
                case 'Ehyy'
                    m_Dat(iEje) = epsilon_Macro(2);
                case 'Ehzz'
                    m_Dat(iEje) = epsilon_Macro(3);
                case 'Ehxy'
                    m_Dat(iEje) = epsilon_Macro(4);
                case 'Shxx'
                    m_Dat(iEje) = sigmaHomog(1);
                case 'Shyy'
                    m_Dat(iEje) = sigmaHomog(2);
                case 'Shzz'
                    m_Dat(iEje) = sigmaHomog(3);
                case 'Shxy'
                    m_Dat(iEje) = sigmaHomog(4);
                case 'Snn'
                    Set = 1;
                    iPG = 6;
                    
                    p_n_tens = e_DatSet(Set).e_DatElem.pointersVAE.p_n_tens;
                    p_condBif= e_DatSet(Set).e_DatElem.pointersVAE.p_condBif ;
                    
                    condBif        =  e_VarAux(Set).VarAuxElem(p_condBif,:) ;
                    %i_stressTilde = e_DatSet(Set).e_DatElem.pointersVHE.i_stressTilde;
                    
                    if condBif
                        
                        index=[(iPG-1)*e_VG.ntens+1:iPG*e_VG.ntens];
                        
                        n_tens= e_VarAux(Set).VarAuxElem(p_n_tens,:);
                        %m_stressTilde=e_VarEst(Set).VarHistElem(i_stressTilde,:);
                        m_stressTilde=e_VarEst(Set).sigma;
                        
                        nx=n_tens(1);  ny=n_tens(4);
                        %tx=-ny; ty=nx;
                        
                        stressTilde=m_stressTilde(index);
                        
                        m_Dat(iEje) = stressTilde(1)*nx^2+2*stressTilde(4)*nx*ny+stressTilde(2)*ny^2 ;
                    else
                        m_Dat(iEje) = 0;
                    end

                case 'Snt'
                    Set = 1;
                    iPG = 6;
                    p_n_tens = e_DatSet(Set).e_DatElem.pointersVAE.p_n_tens;
                    p_condBif= e_DatSet(Set).e_DatElem.pointersVAE.p_condBif ;
                    
                    condBif        =  e_VarAux(Set).VarAuxElem(p_condBif,:) ;
                    %i_stressTilde = e_DatSet(Set).e_DatElem.pointersVHE.i_stressTilde;
                    
                    if condBif
                        
                        index=[(iPG-1)*e_VG.ntens+1:iPG*e_VG.ntens];
                        
                        n_tens= e_VarAux(Set).VarAuxElem(p_n_tens,:);
                        %m_stressTilde=e_VarEst(Set).VarHistElem(i_stressTilde,:);
                        m_stressTilde=e_VarEst(Set).sigma;
                        
                        nx=n_tens(1);  ny=n_tens(4);
                        tx=-ny; ty=nx;
                        
                        stressTilde=m_stressTilde(index);
                        m_Dat(iEje)= stressTilde(1)*nx*tx+stressTilde(4)*(ny*tx+nx*ty)+stressTilde(2)*ny*ty ;
                        
                    else
                        m_Dat(iEje)= 0 ;
                    end
                    
               otherwise
                  error('Archivos de datos: Impresi�n: No est� definido este tipo de dato.')
            end
         end         
%          fId = fopen([fileCompleto,'.cur',num2str(iGraf,'%03d')],'at');
         fId = fopen([fileCompleto,num2str(iGraf,'%03d'),'.dat'],'at');
         fprintf(fId,['%-20.15g',repmat(' %20.15g',1,length(m_Dat)-1),'\n'],m_Dat);
         fclose(fId);            
      end  
   end      
   
end