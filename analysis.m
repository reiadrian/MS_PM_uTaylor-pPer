function analysis(varargin)

% ****************************
% * CONFIGURACION DE ENTORNO *
% ****************************
%Para cerrar el MatLab Pool cuando termina de correr.
cerMPool = 0;

%Definicion de los argumentos de entrada
switch nargin 
   case 0
      %Esto es el caso que no se pasa argumento, donde se se abre el dialogo para abrir el archivo
      %y se utiliza el numero de laboratorios por defecto.
   case 1
      if ischar(varargin{1})
         %Si se pasa como argumento el nombre del archivo con el directorio completo o solo el 
         %nombre si el mismo esta en el path del matlab.
         %analysis(file&Dir)
         [path_file,file,ext] = fileparts(varargin{1});
         file = [file,ext];
      else
         %Si el argumento es un numero, lo toma como el numero de lab a activar en matlabpool.
         %Ver si no verificar el argumento que sea numerico en este caso.
         %analysis(nLabImp)
         nLabImp = varargin{1};         
      end
   case 2
      %Se asume que el primer argumento es un char, indicando el nombre o directorio completo, o
      %solo el directorio.
      if ischar(varargin{2})
         %En el caso que el 2do argumento sea el nombre del archivo.
         %analysis(file,path_file)         
         path_file = varargin{1};
         file = varargin{2};
      else
         %Si el segundo argumento es el numero de laborarios (y el primero es el nombre y directorio
         %del archivo)
         %analysis(file&Dir,nLabImp)
         [path_file,file,ext] = fileparts(varargin{1});
         file = [file,ext];
         nLabImp = varargin{2};   
      end      
   case 3
      %En el caso que se pase los tres argumentos
      %analysis(file,path_file,nLabImp)
      path_file = varargin{1};
      file = varargin{2};
      nLabImp = varargin{3};   
    case 5
      %the case of 4 arguments including the microscale flag
      %analysis(file,path_file,nLabImp)
      path_file = varargin{1};
      file = varargin{2};
      nLabImp = varargin{3};
      isMICRO = varargin{4};
      Snaps = varargin{5};
   otherwise
      error('Analyis: Numero de argumentos de llamada incorrectos.')
end
%En el caso que no se defina el directorio y el nombre del archivo, se busca por medio de un
%dialogo.
if ~exist('file','var') && ~exist('path_file','var')
    [file,path_file] = uigetfile('*.mfl','Definir archivo de calculo');
end    
%Si se quiere imponer una cierta cantidad de laboratorios a crear, se pone un valor distinto mayor que cero.
%Con cero o no se escribe ningun argumento se utiliza el valor por defecto del MatLab o si hay un matlabpool
%abierto, lo utiliza (aunque los labs no sean los por defecto). Un numero negativo hace que no se active el
%matlabpool (si esta abierto, lo cierra) (esto sirve para hacer debug dentro de las funciones que llama el
%parfor, pero no en el mismo parfor, para hacer ello hay que cambiar por el for).
if ~exist('nLabImp','var')
   nLabImp = 0;
end

%Es importante la siguiente linea para que cuando se corta la ejecucion del programa antes de cerrar
%todos los archivos.
fclose all;

p = fileparts(mfilename('fullpath'));
addpath(genpath(p));

% ********************
% * LECTURA DE DATOS *
% ********************
if ischar(file) && ischar(path_file)
   
   fprintf('Inicio de la lectura de datos.\n')
   ticIDLect = tic;
   [in,xx,m_SetElem,f,funbc,e_DatSet,e_VG] = read_data(file,path_file);
   fprintf('Tiempo de lectura y preproceso de datos: %f.\n',toc(ticIDLect));
   
   % *************************************
   % * RESOLUCION DEL PROBLEMA NO LINEAL *
   % *************************************
   %Verificar que la 8.1 es la primera version que aparecio el parcluster.
   if verLessThan('matlab','8.1')
      nLab = matlabpool('size'); %#ok<DPOOL>
      %Se guarda el objeto scheduler para la configuracion activa para ver si la corrida es de tipo local,
      %asi no se realiza la transferencia de archivos.
      objSch = findResource(); %#ok<DFNDR>
      if (nLabImp==0||nLabImp==nLab)&&nLab>0
         disp(['Se utiliza la configuracion activa: ',get(objSch,'configuration'),...
            ' de ',num2str(nLab),' laboratorios.']);
      else
         %Esto tambien cierra el MatLab si nLabImp==-1
         if nLab>0 && nLab~=nLabImp
            matlabpool close %#ok<DPOOL>
         end

         %Recordar que para pocos elementos no conviene paralelizar.
         if nLabImp==0
            if strcmp(objSch.type,'local')
               matlabpool open %#ok<DPOOL>
            else
               matlabpool('open','FileDependencies',{p}) %#ok<DPOOL>
            end
         elseif nLabImp>0
            if strcmp(objSch.type,'local')
               matlabpool('open',nLabImp) %#ok<DPOOL>
            else
               matlabpool('open',nLabImp,'FileDependencies',{p}) %#ok<DPOOL>
            end
         end
      end
      %
      nLab = matlabpool('size'); %#ok<DPOOL>
      if nLab>0 && ~strcmp(objSch.type,'local')
         %Para que funcione en un cluster hay que indicar que archivos se debe copiar (se indica el
         %directorio del programa)
         %matlabpool('addfiledependencies',{p})
         %El directorio por defecto no es el directorio donde copia los archivos, sino que otro dado
         %por pctRunOnAll getFileDependencyDir (directorio donde descomprime los archivos de la
         %dependencia en cada nodo), por lo que en cada nodo se debe modificar para que el directorio
         %activo (pwd) sea el del programa y luego generar los caminos.
         stringCmd = 'pDep=getFileDependencyDir;if ~isempty(pDep),cd(getFileDependencyDir),end;';
         pctRunOnAll(stringCmd)
         %pctRunOnAll pwd
         pctRunOnAll('addpath(genpath(pwd))')
         %pctRunOnAll path
         %Para asegurar que se tenga las ultimas copias del codigo (si no se hace esto cuando se
         %modifica un archivo y el pool se deja abierto, se ignora los cambios, verificar esto)
         matlabpool updatefiledependencies  %#ok<DPOOL>
      end

      %Se almacena los numeros de laborarios y el tipo de conexion para usar imprimir en el newton en un archivo
      %los pasos que no convergen (en el caso que es un cluster tira un error porque no puede crear el archivo en
      %en los nodos, en el camino indicado).
      e_VG.nLab = nLab;
      e_VG.tipoMPool = objSch.type;
      if isfield(objSch,'HasSharedFilesystem')
         e_VG.SharedFS = objSch.HasSharedFilesystem;
      else
         e_VG.SharedFS = false;
      end
   else
      
      o_Pool = [];%gcp('nocreate');
      if isempty(o_Pool)
         nLab = 0;
      else
         nLab = o_Pool.NumWorkers;
      end
      if (nLabImp==0 || nLabImp==nLab) && nLab>0
         disp(['Se utiliza la configuracion activa: ',o_Pool.Cluster.Profile,...
            ' de ',num2str(nLab),' laboratorios.']);
      else
         %Esto tambien cierra el MatLab si nLabImp==-1 o si los Labs a abrir son distintos al que esta.
         if nLab>0 && nLab~=nLabImp
            %Cuidado si esta activado la preferencia que automaticamente cree el pool la opcion -1 no sirve. 
            delete(o_Pool)
         end

         %Recordar que para pocos elementos no conviene paralelizar.
         %o_Cluster = parcluster;
         %Esta version tiene un attach de archivos automatico, que parece que funciona bien, habria que ver si
         %en cada llamada del parfor es mas lento ya que verifica si no se modifico los archivos o los copia
         %de nuevo directamente.
         if nLabImp==0
            %Esto solo sirve para el caso en que el cluster sea MJS.   
            o_Pool = parpool('SpmdEnabled',false);
         elseif nLabImp>0
               o_Pool = parpool(nLabImp,'SpmdEnabled',false);
         end
      end
      %
      if nLabImp<0
         %Arriba ya se asegura que esta cerrado el pool.
         nLab = 0;
      else
         nLab = o_Pool.NumWorkers;
      end
      %
      if nLab>0 && ~strcmp(o_Pool.Cluster.Type,'Local')
         %Para que funcione en un cluster hay que indicar que archivos se debe copiar (se indica el
         %directorio del programa)
         %matlabpool('addfiledependencies',{p})
         %El directorio por defecto no es el directorio donde copia los archivos, sino que otro dado
         %por pctRunOnAll getFileDependencyDir (directorio donde descomprime los archivos de la
         %dependencia en cada nodo), por lo que en cada nodo se debe modificar para que el directorio
         %activo (pwd) sea el del programa y luego generar los caminos.
         %stringCmd = 'pDep=getAttachedFilesFolder;if ~isempty(pDep),cd(getAttachedFilesFolder),end;';
         %pctRunOnAll(stringCmd)
         %pctRunOnAll pwd
         %Parece que la copia automatica de archivos crea directorios separados para cada sub directorio donde
         %esta los archivos .m, y parece que no es necesario agregar al path eso directorios.
         %pctRunOnAll('addpath(genpath(pwd))')
         %pctRunOnAll path
         %Para asegurar que se tenga las ultimas copias del codigo (si no se hace esto cuando se
         %modifica un archivo y el pool se deja abierto, se ignora los cambios, verificar esto)
         o_Pool.updateAttachedFiles
      end

      %Se almacena los numeros de laborarios y el tipo de conexion para usar imprimir en el newton en un archivo
      %los pasos que no convergen (en el caso que es un cluster tira un error porque no puede crear el archivo en
      %en los nodos, en el camino indicado).
      e_VG.nLab = nLab;
      if nLab>0
         e_VG.tipoMPool = o_Pool.Cluster.Type;
          e_VG.SharedFS = false;
      else
         e_VG.tipoMPool = NaN;
         e_VG.SharedFS = false;
      end         
      
   end
   
   e_VG.isMICRO = isMICRO ;
   e_VG.Snap    = Snaps ;
   
   fprintf('Inicio del calculo no lineal.\n')
   ticIDNLA = tic;
   u = nonlinear_analysis(in,xx,m_SetElem,f,funbc,e_DatSet,e_VG); 
   fprintf('**************************************************************\n')
   fprintf('Tiempo de calculo de nonlinear analysis: %f seg.\n',toc(ticIDNLA))
   fprintf('**************************************************************\n')
   
   if cerMPool
      if verLessThan('matlab','8.1')
         matlabpool close
      else
         delete(o_Pool)
      end
   end
   
end

end
