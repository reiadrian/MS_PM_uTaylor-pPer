function e_DatMatME = f_VarME(nomArchME,dirDat,m_ElemPGImpr)

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

   % FUNCI�N CONDICIONES DE BORDE
   %En es caso del modelo multiescala cl�sico se puede mantener las matrices de condiciones de borde
   %en la estructura de propiedades e_DatMat macro ya que no cambian durante el problema.
%    [m_LinCond,vfix,m_InvCRR,doff,dofl,doffCondCte] = f_CondBord(e_VG,xx,...
%      e_DatSet,e_VG.m_ConecFront);

   % Script con variables 
   %(ver si correlo ac�, que significa que hay que llevar varias variables o en cada paso tiempo, donde habr�a
   %que analizar cu�nto tiempo tarda)
   %Ocurre un error al usar usar run, no s� si pasa lo mismo con eval, ya que matlab no se da cuenta que un
   %script o funci�n fue modificada para precompilarla (usa la precompilaci�n anterior). Esto hace que las
   %modificaciones del script las ignora y usa por ejemplo valores de variable que son de la versi�n del
   %script anterior. Esto se arregla con clear all, pero eso puede traer muchos problemas, adem�s que se
   %desaprovecha las precompilaciones previas. Se borra solo la precompilaci�n de la funci�n (hay que usar el
   %nombre de la funci�n solo, sin camino).
   if exist([e_VG.fileCompleto,'.m'],'file')
      clear(nomArchME)
      %run corre el script sin estar en el path o que en el directorio activo
      run([e_VG.fileCompleto,'.m'])
   end
   
   % C�lculo de la medida de la celda unitaria
   %C�mo se hace operaciones sobre el volumen, se guarda el volumen de los elementos en una matriz.
   %Se la ordena seg�n la numeraci�n (orden) de la matriz de conectividad global (como fueron 
   %ingresados los elementos).
   if e_VG.protype==0
                m_VolElem = zeros(1,e_VG.nElem);
                m_VolElem([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem];
   %En caso que se ingrese en el script el �rea micro, se utiliza directamente esa.
                if ~exist('Omega_micro','var')
      %Esta expresi�n es correcta solo si la celda unitaria no tiene agujeros.
                    Omega_micro = sum(m_VolElem);
                end
                e_DatMatME = struct('in',in,'xx',xx,'m_SetElem',m_SetElem,'f',f,'funbc',funbc,...
                'e_DatSet',e_DatSet,'m_ElemPGImpr',m_ElemPGImpr,'omegaMicro',Omega_micro,'e_VG',e_VG);
   elseif e_VG.protype==1
                m_VolElem_d = zeros(1,e_VG.nElem);
                m_VolElem_d([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem_d];
                m_VolElem_p = zeros(1,e_VG.nElem);
                m_VolElem_p([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem_p];
   %En caso que se ingrese en el script el �rea micro, se utiliza directamente esa.
                if ~exist('Omega_micro','var')
      %Esta expresi�n es correcta solo si la celda unitaria no tiene agujeros.
                    Omega_micro_d = sum(m_VolElem_d);
                    Omega_micro_p = sum(m_VolElem_p);
                end
                e_DatMatME = struct('in',in,'xx',xx,'m_SetElem',m_SetElem,'f',f,'funbc',funbc,...
                'e_DatSet',e_DatSet,'m_ElemPGImpr',m_ElemPGImpr,'omegaMicro_d',Omega_micro_d,...
                'omegaMicro_p',Omega_micro_p,'e_VG',e_VG);
   end
   %Imprensi�n de archivo de postprocesado de la malla y inicializaci�n del archivo de datos
   %matlab2gid_mesh(in,xx,e_DatSet,e_VG)
   %f_InicArchDat(e_VG)
     
end