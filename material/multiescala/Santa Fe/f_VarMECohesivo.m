function e_DatMatME = f_VarMECohesivo(nomArchME,dirDat,m_ElemPGImpr,esImplex)

   %Se asume que el archivo de datos de la celda unitaria est� en el mismo directorio el archivo de
   %datos de la malla macro.
   %Los datos de tiempo, ignora los le�dos de este archivo (micro) y usa los le�dos a nivel macro.
   nomArchDat = [nomArchME,'.mfl'];
   [in,xx,m_SetElem,f,funbc,e_DatSet,e_VG] = read_data(nomArchDat,dirDat);
   
   %Para decirle al programa que se est� dentro del c�lculo multiescala
   e_VG.esME = 1;
   %Para impresi�n y debug se guarda la iteraci�n macro (notar que este campo solo se define en el
   %e_VG micro)
   e_VG.iterMacro = [];
   
   % C�lculo de la medida de la celda unitaria
   %C�mo se hace operaciones sobre el volumen, se guarda el volumen de los elementos en una sola matriz.
   %Ver si esto duplica los resultados y mejor obtener la matriz m_VolElem cada vez que se utilice.
   %Se la ordena seg�n la numeraci�n (orden) de la matriz de conectividad global (como fueron 
   %ingresados los elementos).
   m_VolElem = zeros(1,e_VG.nElem);
   m_VolElem([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem];
   %En caso que se ingrese en el script el �rea micro, se utiliza directamente esa.
   if ~exist('Omega_micro','var')
      %Esta expresi�n es correcta solo si la celda unitaria no tiene agujeros.
      Omega_micro = sum(m_VolElem);
   end
   
   % Matriz de elementos vecinos
   m_ElemVec = f_ElemVecino(e_DatSet,e_VG);
   
   e_DatMatME = struct('in',in,'xx',xx,'m_SetElem',m_SetElem,'f',f,'funbc',funbc,'e_DatSet',e_DatSet,...
      'm_ElemPGImpr',m_ElemPGImpr,'m_VolElem',m_VolElem,'omegaMicro',Omega_micro,...
      'm_ElemVec',m_ElemVec,'tipoCBL',[],'pasosPosBif',[],'esImplex',esImplex,'e_VG',e_VG);
   
   % Script con variables 
   %(ver si correlo ac�, que significa que hay que llevar varias variables o en cada paso tiempo,
   %donde habr�a que analizar cu�nto tiempo tarda)
   %Ocurre un error al usar usar run, no s� si pasa lo mismo con eval, ya que matlab no se da cuenta
   %que un script o funci�n fue modificada para precompilarla (usa la precompilaci�n anterior). Esto
   %hace que las modificaciones del script las ignora y usa por ejemplo valores de variable que son
   %de la versi�n del script anterior. Esto se arregla con clear all, pero eso puede traer muchos
   %problemas, adem�s que se desaprovecha las precompilaciones previas. Se borra solo la
   %precompilaci�n de la funci�n (hay que usar el nombre de la funci�n solo, sin camino).
   clear(nomArchME)
   %run corre el script sin estar en el path o que en el directorio activo
   run([e_VG.fileCompleto,'.m'])
   %En este script se define:
   %- m_ElemLoc: indicada como una lista de n�meros de elementos del localizado impuesto.
   %- lMicro: es el ancho de la banda localizada en la celda unitaria. Zona donde se aplica la
   %deformaci�n Bxn/lMicro.
   %- angBanda: el �ngulo de la banda de localizaci�n, que sirve para indentificar cu�l �ngulo se
   %adopta de los dos devueltos del an�lisis de bifurcaci�n.
   %- angNormal: el �ngulo que se quiere forzar que la normal tenga, descartando la que se obtiene
   %del an�lisis de bifurcaci�n.
   %- tipoCBL: el tipo de condici�n de borde del dominio localizado que se impone despu�s de la
   %bifurcaci�n. 0, no se impone CB, 1, se impone m�nima restricci�n, 2, se impone periodicidad (en
   %en este caso habr�a que definir como imponerla, posiblemente crear una funci�n que delvuelva las
   %condiciones de borde), 3, se impone CB lineal, 4, se impone CB de Taylor.
   %- tipoCBBif: se indica el tipo de CB de toda la celda unitaria despu�s de la bifurcaci�n (se
   %debe pasar una funci�n que devuelva m_CondBord redefinido). No implementado todav�a.
   %- pasosPosBif: es la cantidad de pasos despu�s de la bifuraci�n que se quiere imponer el SDA y
   %el modelo multiescala cohesivo.
   %- pasosPosBifCBExt y f_CambioCBExt: Si se quiere cambiar la condici�n externa de la CU despu�s de la
   %bifurcaci�n se debe indicar la cantidad de pasos despu�s de la misma para hacer el cambio y el handle de
   %la funci�n que devuelve la nuevas condiciones de borde.
   %- longFis: longitud de la curva central al dominio localizado (longitud de fisura).
   %- facNormMicro: Media a largo de la fisura S (l�nea central al dominio localizado) de la normales micro.

   if exist('m_ElemLoc','var')
      %Se recupera la matriz de denominaci�n de los elementos (ver si no traerla directamente de
      %read_data).
      m_NumElem = zeros(1,e_VG.nElem);
      for iSet = 1:e_VG.nSet
         m_NumElem(e_DatSet(iSet).m_IndElemSet) = e_DatSet(iSet).m_NumElem;
      end
      %Elementos del dominio localizado impuesto.
      %Es importante poner expl�citamente 1 en any, para el caso que se tenga un solo elemento que
      %localiza (ya que sino el vector lo analiza como una fila).
      e_DatMatME.m_ElemLocImp = any(bsxfun(@eq,m_ElemLoc',m_NumElem),1);
   end
   if exist('lMicro','var')
      %lMicro impuesto
      if isnumeric(lMicro)&&~isnan(lMicro)
         e_DatMatME.lMicroImp = lMicro;
      else
         %Al imponerse lMicro como un string o como NaN se lee las normales micro y las longitudes micro de un
         %archivo de datos.
         e_DatMatME.lMicroImp = NaN;
      end
   end
   if exist('angBanda','var')
      e_DatMatME.angBanda = angBanda;
   end
   if exist('angNormal','var')
      e_DatMatME.angNormalImp = angNormal;
   end
   if exist('tipoCBL','var')
      e_DatMatME.tipoCBL = tipoCBL;
   else
      %Se adopta que por defecto se impone la condici�n de borde lineal en el dominio localizado.
      e_DatMatME.tipoCBL = 3;
   end
   if exist('pasosPosBif','var')
      e_DatMatME.pasosPosBif = pasosPosBif;
   else
      %Por defecto son 0 pasos despu�s de la bifurcaci�n.
      e_DatMatME.pasosPosBif = 0;
   end
   %Para el cambio de la condici�n de borde desp�es de bifurcaci�n
   if exist('pasosPosBifCBExt','var')
      e_DatMatME.pasosPosBifCBExt = pasosPosBifCBExt;
   end
   %Funci�n que devuelve la nueva condici�n de borde externa.
   %En la variable se debe pasar el handle de la funci�n (funciones an�nimas por ejemplo).
   if exist('f_CambioCBExt','var')
      e_DatMatME.f_CambioCBExt = f_CambioCBExt;
   else
      if exist('pasosPosBifCBExt','var')
         error(['Lectura de datos: Propiedades Materiales: Modelo multiescala cohesivo: ',...
            'Si se define la variable pasosPosBifCBExt se debe indicar la funci�n que devuelve las nuevas',...
            'condiciones de borde.'])
      end
   end
   %Longitud total de la fisura en la microescala, se la define si se quiere imponer.
   if exist('longFis','var')
      e_DatMatME.longFisImp = longFis;
   end
   %Media de las normales micro a lo largo de la fisura
   if exist('facNormMicro','var')
      e_DatMatME.facNormMicroImp = facNormMicro;
   end
  
end