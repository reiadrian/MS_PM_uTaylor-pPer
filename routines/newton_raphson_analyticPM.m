function [u,c_GdlCond,Fint2,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT,lambda,o_Par] = newton_raphson_analyticPM(xx,m_LinCond,dofl,doff,u,Du_step_new,...
    c_GdlCond,Du_step_old,Fint2,Fext,e_VarEst_old,e_VarAux,e_DatSet,DefMacro,GradPorMacro,...
    GradPorMacro_up_n,PorMacro,DefMacro_new,GradPorMacro_up_new,PorMacro_new,e_VG,omegaMicro_d)

%******************************************************************************************
%*  ESQUEMA ITERATIVO NEWTON-RAPHSON                                                      *                 
%******************************************************************************************
u_old = u;
% VARIABLES GLOBALES
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
%Se crea aca para evitar estar creandola en cada iteracion del newton. Principalmente para el caso del modelo
%multiescala, donde hay que redifinir hvar como una estructura y no como una matriz como los otros casos. No
%interesa sus valores, solo que esta definida para indexar correctamente los argumentos de salida.
%Ver si no cambiar la estructura para que se pueda indicar que se necesita como variable de estado (sigma,
%eps, eps_fluct, etc.) por set. Esto evitaria tener matrices que no se necesita en ciertos sets. Otra es tener
%variables que indiquen las dimensiones que se necesita de estas matrices, y en el caso que no se necesite que
%sea cero.
%###################################################################################
%Agregue Fext
e_VarEst_new = f_eVarEstInic({'eps','phi','porpr','sigmaE','sigmaT','mflu','velflu',...
        'velflu_sta','velflu_total','eps_fluct','phi_fluct','p_fluct','hvar','bodyForce'},e_DatSet,e_VG,0);
%###################################################################################
du_iter = zeros(ndoft,1);

% ESQUEMA ITERATIVO DE NEWTON
% Variables condensadas
c_GdlCond = f_cVar_initIncNR(e_DatSet,e_VG,c_GdlCond);

% ITERACION "0", EVALUACION DE LA FUERZA INTERNA
e_VG.iter = iter;
%Para forzar la iteracion cero, se resuelva como un problema elastico, independiente de
%como venga la variable e_VG.elast (luego en las siguientes iteraciones se recupera el valor 
%impuesto en esa variable). Esta resuelve el caso que durante la parte elastica del problema, se
%itere dos veces ya que el valor inicial fuerza un comportamiento inelastico ficticio.
%Notar que lo se consigue es realizar en la interacion 1 un avance elastico, es decir 
%u(1) = u(0)+Du (u(i+1)=u(i)+Du), Du elastico (cualquiera sea el estado actual del problema,
%que no corresponde al comportamiento de los modelos constitutivos elegidos) y que el residuo va
%corresponder a ese comportamiento elastico asumido (es decir va ser incorrecto).

%CALCULO DE bj
f_sin = @(y,j) sin(2*j*pi
bj = GradPorMacro_up_n*
%CALULO DE LAS MICRO-PORO PRESIONES FLUCTUANTES


[KT,Fint,u,c_GdlCond,e_VarEst_new,e_VarAux] = f_Par_MatGlobales_analyticPM(xx,u,u_old,c_GdlCond,e_DatSet,...
    e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,GradPorMacro,GradPorMacro_up_n,PorMacro,...
    DefMacro_new,GradPorMacro_up_new,PorMacro_new,e_VG);

% CALCULO DEL RESIDUO (se considera los residuos de los gdl sin condensar)
res = residuo(Fext,Fint,m_LinCond,dofl,doff,lambda);

% NORMA DEL RESIDUO (se considera los residuos de los gdl sin condensar)
norm_res = norm(res)+f_cVar_normRes(e_DatSet,e_VG,c_GdlCond); 

% IMPRESION de interacion 0
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
norm_desp = 1;
norm_por =1;
tol_rel = tolnr;
max_iter = iterMax;

% LOOP DE NEWTON
%Se usa < y no <= en iter < max_iter porque el valor inicial de iter es 0.
while (norm_res>tol_rel||norm_desp>tol_rel||norm_por>tol_rel)&&iter<max_iter||iter<1
   
   % CONTEO DE ITERACIONES
   iter = iter+1;
   e_VG.iter = iter;
   
   % INCREMENTO DE DESPLAZAMIENTO S/ESTRATEGIA DE CONTROL

   [du_iterl,dlambda] = incremental_disp(res,KT,Fext,m_LinCond,dofl,doff,Du_step_new,Du_step_old,e_VG);
   du_iter(dofl) = du_iterl;
   du_iter(doff) = m_LinCond*du_iterl;
   
   % INCREMENTO DE VARIABLES CONDENSADAS (gdl internos)
   c_GdlCond = f_cVar_IncrNR(du_iter,c_GdlCond,e_DatSet,e_VG);
   
   % DESPLAZAMIENTO TOTAL
   u = u + du_iter;
   Du_step_new = Du_step_new + du_iter;
   
   lambda = lambda + dlambda;
   
   if e_VG.CONTROL_STRAT==4
       u(doff) =  e_VG.vfix_doff*lambda;
   end 
   
   % MATRIZ DE RIGIDEZ y COMPUTO DE LA FUERZA INTERNA
   [KT,Fint,u,c_GdlCond,e_VarEst_new,e_VarAux,c_CT,o_Par] =...
        f_Par_MatGlobales_analyticPM(xx,u,u_old,c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,...
        e_VarAux,DefMacro,GradPorMacro,GradPorMacro_up_n,PorMacro,...
        DefMacro_new,GradPorMacro_up_new,PorMacro_new,e_VG);
   
   
%    if protype==1 %AA
%        gdl_carga=e_VG.gdl_carga;
%        Reac=-KT*u;
%        Reac(gdl_carga)=0.0;
%        F=Fext-Reac;
%        r=F-Fint;
%    end

    res = residuo(Fext,Fint,m_LinCond,dofl,doff,lambda);

   
   % NORMA DEL RESIDUO (se considera los residuos de los gdl sin condensar)
   norm_res = norm(res)+f_cVar_normRes(e_DatSet,e_VG,c_GdlCond);
   % NORMA DEL DESPLAZAMIENTO (se considera los residuos de los gdl sin condensar)
   v_desp = (1:ndoft)';
   v_desp(3:3:end) = 0;
   du_despl = du_iter(v_desp.*dofl > 0);
   norm_desp = norm(du_despl);
   % NORMA DE LAS POROPRESIONES (se considera los residuos de los gdl sin condensar)
   v_por = (1:ndoft)';
   v_por(1:3:end) = 0;
   v_por(2:3:end) = 0;
   du_por = du_iter(v_por.*dofl > 0);
   du_por = du_por(du_por ~= 0);
   norm_por = norm(du_por);
   
   % IMPRESION POR PANTALLA DEL PROCESO DE CONVERGENCIA
   if ~e_VG.esME
       disp('=================================================================');
       fprintf('ITERACION: %2d (ITERACION MACRO)\n',iter);
       fprintf('|res| = %-8.3e (ITERACION MACRO)\n',norm_res);
       fprintf('|Ddesp| = %-8.3e (ITERACION MACRO)\n',norm_desp);
       fprintf('|Dpore| = %-8.3e (ITERACION MACRO)\n',norm_por);
       disp('=================================================================');
   else
       fprintf('ITERACION: %2d ; \n',iter);
       fprintf('|res| = %-8.3e \n',norm_res);
       fprintf('|Ddesp| = %-8.3e \n',norm_desp);
       fprintf('|Dpore| = %-8.3e \n',norm_por);
   end
end

Fint2 = Fint2 + Fint;

% IMPRESION POR PANTALLA DEL PROCESO DE CONVERGENCIA
if iter==max_iter&&norm_res>tol_rel
   %En el debug se puede hacer que se active un breakpoint si usa un warning como mensaje.
   warning('off','backtrace')
   if e_VG.esME
      %Para tener en cuenta que si esta paralelizado la macroescala y no tiene carpeta compartidas, el newton
      %de la microescala no va encontrar donde guardar este archivo, por lo que en ese caso se desactiva.
      if e_VG.nLab==0||strcmp(e_VG.tipoMPool,'local')||e_VG.SharedFS
         nomArch = [e_VG.fileCompleto,'.pasNoConv'];
         fId = fopen(nomArch,'at');
         fprintf(fId,'Paso %d - IterMacro %d - ElemMacro %d - PGMacro %d - Iter %d: |res| = %e\n',e_VG.istep,...
            e_VG.iterMacro,e_VG.iElemNumMacro,e_VG.iPGMacro,iter,norm_res);
         fclose(fId);
      end
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
end