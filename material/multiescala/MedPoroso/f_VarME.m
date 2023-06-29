function e_DatMatME = f_VarME(nomArchME,dirDat,m_ElemPGImpr)

   %Se asume que el archivo de datos de la celda unitaria esta en el mismo directorio el archivo de
   %datos de la malla macro.
   %Los datos de tiempo, ignora los leidos de este archivo (micro) y usa los leidos a nivel macro.
   nomArchDat = [nomArchME,'.mfl'];
   [in,xx,m_SetElem,f,funbc,e_DatSet,e_VG] = read_data(nomArchDat,dirDat);
   
   %Para decirle al programa que se esta dentro del calculo multiescala
   e_VG.esME = 1;

   %Para impresion y debug se guarda la iteracion macro (notar que este campo solo se define en el
   %e_VG micro)
   e_VG.iterMacro = [];

   % FUNCION CONDICIONES DE BORDE
   %En es caso del modelo multiescala clasico se puede mantener las matrices de condiciones de borde
   %en la estructura de propiedades e_DatMat macro ya que no cambian durante el problema.
%    [m_LinCond,vfix,m_InvCRR,doff,dofl,doffCondCte] = f_CondBord(e_VG,xx,...
%      e_DatSet,e_VG.m_ConecFront);

   % Script con variables 
   %(ver si correlo aca, que significa que hay que llevar varias variables o en cada paso tiempo, donde habria
   %que analizar cuanto tiempo tarda)
   %Ocurre un error al usar usar run, no si si pasa lo mismo con eval, ya que matlab no se da cuenta que un
   %script o FUNCION fue modificada para precompilarla (usa la precompilacion anterior). Esto hace que las
   %modificaciones del script las ignora y usa por ejemplo valores de variable que son de la version del
   %script anterior. Esto se arregla con clear all, pero eso puede traer muchos problemas, ademas que se
   %desaprovecha las precompilaciones previas. Se borra solo la precompilacion de la FUNCION (hay que usar el
   %nombre de la FUNCION solo, sin camino).
   if exist([e_VG.fileCompleto,'.m'],'file')
      clear(nomArchME)
      %run corre el script sin estar en el path o que en el directorio activo
      run([e_VG.fileCompleto,'.m'])
   end
   
   % calculo de la medida de la celda unitaria
   %Como se hace operaciones sobre el volumen, se guarda el volumen de los elementos en una matriz.
   %Se la ordena segun la numeracion (orden) de la matriz de conectividad global (como fueron 
   %ingresados los elementos).
   if e_VG.protype==0
                m_VolElem = zeros(1,e_VG.nElem);
                m_VolElem([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem];
   %En caso que se ingrese en el script el area micro, se utiliza directamente esa.
                if ~exist('Omega_micro','var')
      %Esta expresion es correcta solo si la celda unitaria no tiene agujeros.
                    Omega_micro = sum(m_VolElem);
                end
                e_DatMatME = struct('in',in,'xx',xx,'m_SetElem',m_SetElem,'f',f,'funbc',funbc,...
                'e_DatSet',e_DatSet,'m_ElemPGImpr',m_ElemPGImpr,'omegaMicro',Omega_micro,'e_VG',e_VG);
    %########################################################################################
   elseif e_VG.protype==1
                m_VolElem_d = zeros(1,e_VG.nElem);
                m_VolElem_d([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem_d];
                m_VolElem_p = zeros(1,e_VG.nElem);
                m_VolElem_p([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem_p];
   %En caso que se ingrese en el script el area micro, se utiliza directamente esa.
                if ~exist('Omega_micro','var')
      %Esta expresion es correcta solo si la celda unitaria no tiene agujeros.
                    Omega_micro_d = sum(m_VolElem_d);
                    Omega_micro_p = sum(m_VolElem_p);
                end
                e_DatMatME = struct('in',in,'xx',xx,'m_SetElem',m_SetElem,'f',f,'funbc',funbc,...
                'e_DatSet',e_DatSet,'m_ElemPGImpr',m_ElemPGImpr,'omegaMicro_d',Omega_micro_d,...
                'omegaMicro_p',Omega_micro_p,'e_VG',e_VG);
    %########################################################################################
   end
   %Imprension de archivo de postprocesado de la malla y inicializacion del archivo de datos
   %matlab2gid_mesh(in,xx,e_DatSet,e_VG)
   %f_InicArchDat(e_VG)
     
end