function [u,c_GdlCond,Fint2,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT,lambda,o_Par] = newton_raphson(xx,...
   m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint2,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
   DefMacro,e_VG)

%******************************************************************************************
%*  ESQUEMA ITERATIVO NEWTON-RAPHSON                                                      *
%*                                                                                        *
%*  ARGUMENTOS DE ENTRADA:                                                                *
%*  inn         : lista de nodos                                                          *
%*  xx          : lista de coordenadas                                                    *
%*  conec       : lista de conectividades                                                 *
%*  u           : vector de desplazamientos                                               *                  
%*  Du_step_old : incremento de desplazamiento en el paso previo                          *
%*  Fint        : vector de fzas internas                                                 *
%*  Fext        : vector de fzas externa                                                  *
%*  hvar_old    : arreglo de variables historicas del paso previo                         *
%*  Eprop       : lista de propiedades de los elementos                                   *
%*                                                                                        *                  
%*  ARGUMENTOS DE SALIDA:                                                                 *                  
%*  u           : vector de desplazamientos                                               *                  
%*  Fint        : vector de fzas internas                                                 *
%*  sigma_new   : arreglo de tensiones del paso actual                                    *
%*  epsr_new    : arreglo de deformaciones del paso actual                                *
%*  hvar_new    : arreglo de variables historicas del paso actual                         *
%*  Du_step_new : incremento de desplazamiento en el paso actual                          *
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% VARIABLES GLOBALES
protype = e_VG.protype; %AA
ndoft = e_VG.ndoft;
tolnr = e_VG.tolnr;
iterMax = e_VG.iterMax;
if e_VG.CONTROL_STRAT==4
 lambda = e_VG.lambda;
else
 lambda=1;
end

% INICIALIZACION DE VARIABLES
iter = 0;
%Se crea ac� para evitar estar cre�ndola en cada iteraci�n del newton. Principalmente para el caso del modelo
%multiescala, donde hay que redifinir hvar como una estructura y no como una matriz como los otros casos. No
%interesa sus valores, solo que est� definida para indexar correctamente los argumentos de salida.
%Ver si no cambiar la estructura para que se pueda indicar que se necesita como variable de estado (sigma,
%eps, eps_fluct, etc.) por set. Esto evitar�a tener matrices que no se necesita en ciertos sets. Otra es tener
%variables que indiquen las dimensiones que se necesita de estas matrices, y en el caso que no se necesite que
%sea cero.
if protype==0 %AA
e_VarEst_new = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','bodyForce'},e_DatSet,e_VG,0);
elseif protype==1 %AA
e_VarEst_new = f_eVarEstInic({'sigmaE','sigmaT','eps','eps_fluct','hvar','bodyForce'},e_DatSet,e_VG,0);
end %AA

du_iter = zeros(ndoft,1);

% ESQUEMA ITERATIVO DE NEWTON
% Variables condensadas
c_GdlCond = f_cVar_initIncNR(e_DatSet,e_VG,c_GdlCond);

% MODIFICACI�N DEL VALOR INICIAL
%Parece que es muy lento, considerando que en multiescala se lo llama repetida veces (una vez por
%cada iteraci�n macro). Ver si se calcula esta condici�n inicial, no utilizar la iteraci�n cero y
%forzar que siempre se haga una iteraci�n m�s.
%ticIdVI = tic;
%u = f_ValorInicial(xx,m_LinCond,dofl,doff,u,c_GdlCond,Fext,e_VarEst_new,e_VarEst_old,...
%   e_VarAux,e_DatSet,DefMacro,e_VG);
%fprintf('Tiempo de c�lculo de valor inicial: %f\n',toc(ticIdVI));

% ITERACION "0", EVALUACION DE LA FUERZA INTERNA
e_VG.iter = iter;
%Para forzar la iteraci�n cero, se resuelva como un problema el�stico, independiente de
%como venga la variable e_VG.elast (luego en las siguientes iteraciones se recupera el valor 
%impuesto en esa variable). Esta resuelve el caso que durante la parte el�stica del problema, se
%itere dos veces ya que el valor inicial fuerza un comportamiento inel�stico ficticio.
%Notar que lo se consigue es realizar en la interaci�n 1 una avance el�stico, es decir 
%u(1) = u(0)+Du (u(i+1)=u(i)+Du), Du el�stico (cualquiera sea el estado actual del problema,
%que no corresponde al comportamiento de los modelos constitutivos elegidos) y que el residuo va
%corresponder a ese comportamiento el�stico asumido (es decir va ser incorrecto).
% elastOrig = e_VG.elast;
% e_VG.elast = 1;
%
%ticIDEnsam = tic;
[KT,Fint,c_GdlCond,e_VarEst_new,e_VarAux] = f_Par_MatGlobales(xx,u,Du_step_new,c_GdlCond,e_DatSet,...
   e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,e_VG); %AA
%fprintf('Tiempo de ensamblaje Iteraci�n 0: %f\n',toc(ticIDEnsam));
% e_VG.elast = elastOrig;
Fint2 = Fint2 + Fint;
% ITERACION "0", COMPUTO DEL RESIDUO
res = residuo(Fext,Fint2,m_LinCond,dofl,doff,lambda);
% NORMA DEL RESIDUO (se considera los residuos de los gdl sin condensar)
norm_res = norm(res)+f_cVar_normRes(e_DatSet,e_VG,c_GdlCond); %AA: add este desde abajo
%norm_res =
%norm([res;cell2mat((cellfun(@(x)reshape(x,[],1),c_GdlCond(:,2),'UniformOutPut',0)))],inf);%Bloque
%este

% IMPRESION de interaci�n 0
if ~e_VG.esME
   disp('=================================================================');
   fprintf('ITERACION: %2d ; |res| = %-8.3e (VALOR INICIAL-ITERACION MACRO)\n',iter,norm_res);
   disp('=================================================================');
else
   fprintf('Problema micro: Elemento Macro %d: Punto gauss macro %d:\n',e_VG.iElemNumMacro,e_VG.iPGMacro);
   fprintf('ITERACION: %2d ; |res| = %-8.3e (VALOR INICIAL)\n',iter,norm_res);
end

% INICIALIZACION DE VARIABLES PROPIAS DE ESQUEMA NEWTON-RAPHSON
norm_res = 1;
tol_rel = tolnr;
max_iter = iterMax;

% LOOP DE NEWTON
%Se usa < y no <= en iter < max_iter porque el valor inicial de iter es 0.
while norm_res>tol_rel&&iter<max_iter||iter<1
   
   % CONTEO DE ITERACIONES
   iter = iter+1;
   e_VG.iter = iter;
   
   % INCREMENTO DE DESPLAZAMIENTO S/ESTRATEGIA DE CONTROL
   %ticIDResSist = tic;
   [du_iterl,dlambda] = incremental_disp(res,KT,Fext,m_LinCond,dofl,doff,Du_step_new,Du_step_old,e_VG);
   du_iter(dofl) = du_iterl;
   du_iter(doff) = m_LinCond*du_iterl;
   %fprintf('Tiempo de resoluci�n del sistema: %f\n',toc(ticIDResSist));
   
   % INCREMENTO DE VARIABLES CONDENSADAS (gdl internos)
   c_GdlCond = f_cVar_IncrNR(du_iter,c_GdlCond,e_DatSet,e_VG);
   
   % DESPLAZAMIENTO TOTAL
   u = u + du_iter;
   Du_step_new = Du_step_new + du_iter;
   
   lambda = lambda + dlambda;
   
   if e_VG.CONTROL_STRAT==4
       u(doff) =  e_VG.vfix_doff*lambda;
   end 

   % VARIABLES CONDESADAS TOTALES
 %  c_GdlCond(:,1) = cellfun(@plus,c_GdlCond(:,1),c_dGdlCond,'UniformOutput',false);
   
   % MATRIZ DE RIGIDEZ y COMPUTO DE LA FUERZA INTERNA
   %ticIDEnsam = tic;
   [KT,Fint,c_GdlCond,e_VarEst_new,e_VarAux,c_CT,o_Par] = f_Par_MatGlobales(xx,u,Du_step_new,c_GdlCond,...
      e_DatSet,e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,e_VG);
   %fprintf('Tiempo de ensamblaje: %f\n',toc(ticIDEnsam));
   Fint2 = Fint2+Fint;
%    if protype==1 %AA
%        gdl_carga=e_VG.gdl_carga;
%        Reac=-KT*u;
%        Reac(gdl_carga)=0.0;
%        F=Fext-Reac;
%        r=F-Fint;
%    end
   
%    if iter==1 && e_VG.istep==1 %AA 
%        Fint=0.95*Fint; %Obligo a no converger en la primera iteracion
%    end %AA
   
   % EVALUACION DEL RESIDUO
   res = residuo(Fext,Fint2,m_LinCond,dofl,doff,lambda);
        
   % NORMA DEL RESIDUO (se considera los residuos de los gdl sin condensar)
   norm_res = norm(res)+f_cVar_normRes(e_DatSet,e_VG,c_GdlCond);
   %norm_res = norm([res;cell2mat((cellfun(@(x)reshape(x,[],1),c_GdlCond(:,2),...
   %  'UniformOutPut',0)))],inf);
   
   % IMPRESION POR PANTALLA DEL PROCESO DE CONVERGENCIA
    if ~e_VG.esME
        disp('=================================================================');
        fprintf('ITERACION: %2d ; |res| = %-8.3e (ITERACION MACRO)\n',iter,norm_res);
        disp('=================================================================');
    else
        fprintf('ITERACION: %2d ; |res| = %-8.3e\n',iter,norm_res);
    end
    
end

% IMPRESION POR PANTALLA DEL PROCESO DE CONVERGENCIA
if iter==max_iter&&norm_res>tol_rel
   %fprintf('ESQUEMA DE NEWTON NO CONVERGE\n');
   %pause
   %En el debug se puede hacer que se active un breakpoint si usa un warning como mensaje.
   warning('off','backtrace')
   if e_VG.esME
      %Para tener en cuenta que si est� paralelizado la macroescala y no tiene carpeta compartidas, el newton
      %de la microescala no va encontrar donde guardar este archivo, por lo que en ese caso se desactiva.
      if e_VG.nLab==0||strcmp(e_VG.tipoMPool,'local')||e_VG.SharedFS
         nomArch = [e_VG.fileCompleto,'.pasNoConv'];
         fId = fopen(nomArch,'at');
         fprintf(fId,'Paso %d - IterMacro %d - ElemMacro %d - PGMacro %d - Iter %d: |res| = %e\n',e_VG.istep,...
            e_VG.iterMacro,e_VG.iElemNumMacro,e_VG.iPGMacro,iter,norm_res);
         fclose(fId);
      end
      %
      warning(['Problema micro: Elemento Macro %d: Punto de gauss macro %d: Esquema de newton no converge ',...
         '(|res|=%e).'],e_VG.iElemNumMacro,e_VG.iPGMacro,norm_res);
   else
      nomArch = [e_VG.fileCompleto,'.pasNoConv'];
      fId = fopen(nomArch,'at');
      fprintf(fId,'Paso %d - Iter %d: |res| = %e\n',e_VG.istep,iter,norm_res);
      fclose(fId);
      warning('ESQUEMA DE NEWTON NO CONVERGE'); %#ok<WNTAG>
   end
   warning('on','backtrace')
else
   if ~e_VG.esME
      fprintf('CONVERGE\n');
   else
      fprintf('Problema micro: Converge\n')
   end
   %fprintf('----------------------------------\n');
end