function u = nonlinear_analysis(in,xx,m_SetElem,f,funbc,e_DatSet,e_VG)

%***************************************,**************************************************
%*  RUTINA PARA RESOLVER EL PROBLEMA NO LINEAL MATERIAL MEDIANTE M.E.F.                   *
%******************************************************************************************

% VARIABLES GLOBALES
struhyp = e_VG.struhyp;
ndoft = e_VG.ndoft;
nElem = e_VG.nElem;
ndn = e_VG.ndn;
np = e_VG.np;
postpro_impre_step = e_VG.postpro_impre_step;
CONTROL_STRAT = e_VG.CONTROL_STRAT;
ndime = e_VG.ndime;
ntens = e_VG.ntens;
nSet = e_VG.nSet;

if e_VG.protype==0 %AA: Dtime variable en problemas bifase. P/monofase es constante 
    Dtime = e_VG.Dtime;
end

% INICIALIZACION DE VARIABLES

Fint = zeros(ndoft,1);          % Vector de fuerzas internas
Fext = zeros(ndoft,1);          % Vector de fuerzas externas

Fint2 = zeros(ndoft,1);          % Vector de fuerzas internas
Fext2 = zeros(ndoft,1);          % Vector de fuerzas externas
u = zeros(ndoft,1);             % Vector de desplazamientos totales
Du_step_old = zeros(ndoft,1);   % Vector de incrementos de desplazamientos en el paso de tiempo previo
%Por simplicidad se considera una deformacion macro aplicada distinta por elemento, y no por PG.
%(no es necesario considerar una estructura para esta variable)
%DefMacro = zeros(ntens*npg,nElem);
%DefMacro = zeros(ntens,nElem);
sigmaHomog = zeros(ntens,1); % AA: no deberia ser (ntens,numero elemntos micro/macro)

if e_VG.isMICRO.MICRO
    %Para analisis RVE
    switch struhyp %e_VG.conshyp
        case 20 %{100,110} % No linealidad geometrica Elastic Material neo-Hookean y J2, It is assumed next order Voigt tensor [Fxx;Fyy;Fzz;Fxy;Fyx]
            epsilon_Macro  = [1; 1; 1; 0; 0;];
            epsilon_Macro0 = e_VG.isMICRO.epsilon_Macro0 ;% +epsilon_Macro ;% % Increment of macro-strain
        otherwise % Pequeñas deformaciones
            epsilon_Macro  = zeros(ntens,1); % Initial imposed macro-strain
            epsilon_Macro0 = e_VG.isMICRO.epsilon_Macro0 ; % Increment of macro-strain
            %Caso bifasico (suelo saturado)
            if e_VG.protype==1
                p_Macro  = zeros(1,1); %Macro-poro presiones iniciales impuestas
                p_Macro0 = e_VG.isMICRO.p_Macro0; %Incremento de macro-poro presiones
                phi_Macro  = zeros(2,1); %Macro-gradiente de poro-presiones iniciales impuestas
                phi_Macro0= e_VG.isMICRO.phi_Macro0; %Mac                
            end
    end
else
    %Para analisis MULTIESCALA
    epsilon_Macro = zeros(ntens,1); %JLM %Variable necesaria en MULTIESCALA solo para impresion de resultados
end

%!!!! Iniciacion Matriz de Snapshots (JLM)
if e_VG.isMICRO.MICRO && ~isempty(e_VG.Snap)
    [dir,nomb] = fileparts(e_VG.fileCompleto);
    nglT = 0 ;
    ntens= e_VG.ntens ;  
    for iset=1:nSet
        nElemS = e_DatSet(iset).nElem;
        npg    = e_DatSet(iset).e_DatElem.npg;
        nglT   = nglT + nElemS*npg;
    end
    PointersToSet1 = true(nglT*ntens,1);   % Elementos del dominio regular
    PointersToSet2 = false(nglT*ntens,1);  % Elementos del dominio singular
    iSetDis = 1;
    if nSet>1
        iSetDis = 2;
    end
    
    % SOLO FUNCIONA CON 2 subdominios
    ind_elem_set = [];
    nElemS = e_DatSet(iSetDis).nElem; 
    npg = e_DatSet(iSetDis).e_DatElem.npg;
    ElemSet = e_DatSet(iSetDis).m_IndElemSet;

    for ielem=1:nElemS
        ind_elem_set = [ind_elem_set ; [(ElemSet(ielem)-1)*npg*ntens + [1:1:npg*ntens]]'];   
    end
    PointersToSet1(ind_elem_set,1) = false;
    PointersToSet2(ind_elem_set,1) = true; 
    save([dir '/DomainPointers.mat'],'PointersToSet1','PointersToSet2')
    clear ElemSet PointersToSet1 PointersToSet1
    %     Estructura de datos para almacenar los snapshots
    Snapshots = struct('SnapStrain',[],'SnapStress',[],'SnapWeight',[],...
        'SnapEnergy_e',[],'SnapEnergy_e_vol',[],'SnapEnergy_e_dev',[],...
        'SnapEnergy_p',[],'SnapEnergy_t',[],'Snapflag',[]);
end

ELOCCalc = false(1,nElem);
if e_VG.protype==0
    e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar'},e_DatSet,e_VG); %AA: add case 14
elseif e_VG.protype==1 
%####################################################################################
    e_VarEst_old = f_eVarEstInic({'eps','phi','porpr','sigmaE','sigmaT','mflu','velflu',...
        'velflu_sta','velflu_total','eps_fluct','phi_fluct','p_fluct','hvar','bodyForce'},e_DatSet,e_VG);
%####################################################################################
end

%Se envia xx solo para inicializar los matriz de deformacion del punto central del elemento FBar_LD.
%####################################################################################
e_VarAux = f_eVarAuxInic(xx,e_DatSet,e_VG);
%####################################################################################
c_GdlCond = f_cVarCondInic(e_DatSet,e_VG);

%Impresion de archivo de postprocesado de la malla e inicializacion del archivo de datos
matlab2gid_mesh(in,xx,e_DatSet,e_VG)
e_VG.istep = 0;
f_InicArchDat(in,m_SetElem,e_DatSet,e_VG)

f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_old,e_VarAux,e_DatSet,m_SetElem,sigmaHomog,epsilon_Macro,e_VG) ; 

% VALOR DE LA FUNCION Psi EN EL TIEMPO CERO
psi_value_old      = get_psi_value(funbc,0,e_VG);

%Con istepSave se indica desde que paso se empieza a correr. Si lee el archivo mat, el valor nulo de
%istepSave es pisado por el paso donde se guarda el workspace.
istepSave = 0;
[dirFile,nomFile] = fileparts(e_VG.fileCompleto);

%Se indica cada cuanto se salva los pasos
deltaPasoSave = 100;

% INFORMACION ADICIONAL
%Se lee un script de matlab. Esto permite cambiar algunas propiedades, matrices, etc., antes de entrar en el
%calculo en forma rapida sin modificar demasiado el programa.
%Ocurre un error al usar usar run, no se si pasa lo mismo con eval, ya que matlab no se da cuenta
%que un script o funcion fue modificada para precompilarla (usa la precompilacion anterior). Esto
%hace que las modificaciones del script las ignora y usa por ejemplo valores de variable que son de
%la version del script anterior.
%Esto se arregla con clear all, pero eso puede traer muchos problemas, ademas que se desaprovecha
%las precompilaciones previas. Se borra solo la precompilacion de la funcion (hay que usar el nombre
%de la funcion solo, sin camino).
if exist([e_VG.fileCompleto,'.m'],'file')
   clear(nomFile)
   %run corre el script sin estar en el path o que en el directorio activo
   run([e_VG.fileCompleto,'.m'])
end

% FUNCION CONDICIONES DE BORDE
%####################################################################################
[m_LinCond,vfix,m_InvCRR,doff,dofl,doffCondCte] = f_CondBord(e_VG,xx,e_DatSet,e_VG.m_ConecFront);
%####################################################################################

% Medicion de tiempos en un cluster mediante parTicToc
c_ParTT = cell(np,1);
e_VG.Dtime2 = 0;
%$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#
% e_VG.time_n = 0;
%$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#
% INTEGRACION TEMPORAL
for istep = istepSave+1:np
   
   ticStep = tic;
   
   e_VG.istep = istep;
   %####################################################################################
   if e_VG.protype==0 %AA: Dtime variable en problemas bifase
       time = Dtime*istep;
   elseif e_VG.protype==1 %AA: Dtime variable en problemas bifase
       time=funbc(istep,2);
       e_VG.Dtime = time;
       e_VG.Dtime2 = e_VG.Dtime2 + time;
   end %AA
   %####################################################################################
   disp('=================================================================');
   fprintf('STEP: %-3d\n',istep);
   disp('=================================================================');
   
   % VALOR ACTUAL DE LA FUNCION TEMPORAL PARA ESCALAR CONDICIONES DE BORDE
   psi_value = get_psi_value(funbc,time,e_VG);
   
   %
   if e_VG.isMICRO.MICRO 
       %Para analisis RVE
       epsilon_Macro = epsilon_Macro + epsilon_Macro0*(psi_value - psi_value_old);
       %if e_VG.MOD_TYPE==2 % MONOSCALE + HROM
%       DefMacro =arrayfun(@(x)repmat(epsilon_Macro,[1,x.e_DatElem.npg,x.nElem]),e_VG.ROM_II.e_DatSet_ROM,'UniformOutput',false);
%       else % MONOSCALE + ROM & FE
       DefMacro =arrayfun(@(x)repmat(epsilon_Macro,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);  
       
       %end
       DefMacroT = epsilon_Macro;  
       %Caso bifasico (suelo saturado)
       if e_VG.protype==1
            p_Macro = p_Macro + p_Macro0*(psi_value - psi_value_old);
            phi_Macro = phi_Macro + phi_Macro0*(psi_value - psi_value_old);
            
            GradPorMacro = arrayfun(@(x)repmat(phi_Macro,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
            GradPorMacro_up_n = arrayfun(@(x)repmat(phi_Macro,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);  
            PorMacro = arrayfun(@(x)repmat(p_Macro,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
            
            GradPorMacroT = zeros(2,1); %AA22
            PorMacroT = zeros(1,1); %AA22   
       end
   %     
   else
       %Para analisis MULTIESCALA
       switch struhyp %e_VG.conshyp
            case 20 %{100,110} % No linealidad geometrica Elastic Material neo-Hookean y J2
                %It is assumed next order Voigt tensor [Fxx;Fyy;Fzz;Fxy;Fyx]
                DefMacroT = zeros(ntens,nElem);
                DefMacro = arrayfun(@(x)[ones(3,x.nElem);zeros(ntens-3,x.nElem)],e_DatSet,'UniformOutput',false);
%                     DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.e_DatElem.npg,x.nElem),e_DatSet,'UniformOutput',false);
%                     DefMacroT = zeros(ntens,1);
            otherwise % Pequeñsas deformaciones
                %DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.e_DatElem.npg,x.nElem),e_DatSet,'UniformOutput',false);
                DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.nElem),e_DatSet,'UniformOutput',false); %AA
                DefMacro_new = arrayfun(@(x)zeros(e_VG.ntens,x.nElem),e_DatSet,'UniformOutput',false); %AA
                 %AA: Quite la cantidad de PG suponiendo que la deformacion macro es constante en todo el elento y aplicado a cada PG
                DefMacroT = zeros(ntens,1);
                %Caso bifasico (suelo saturado)
                if e_VG.protype==1 
                   GradPorMacro = arrayfun(@(x)zeros(2,x.nElem),e_DatSet,'UniformOutput',false); %AA22 VER DIMENSIONES DEL VECTOR
                   GradPorMacro_up_n = arrayfun(@(x)zeros(2,x.nElem),e_DatSet,'UniformOutput',false); %AA22 VER DIMENSIONES DEL VECTOR  
                   PorMacro = arrayfun(@(x)zeros(1,x.nElem),e_DatSet,'UniformOutput',false); %AA22 VER DIMENSIONES DEL VECTOR
                   GradPorMacro_up_new = arrayfun(@(x)zeros(2,x.nElem),e_DatSet,'UniformOutput',false); %AA22 VER DIMENSIONES DEL VECTOR  
                   PorMacro_new = arrayfun(@(x)zeros(1,x.nElem),e_DatSet,'UniformOutput',false); %AA22 VER DIMENSIONES DEL VECTOR
                   GradPorMacroT = zeros(2,1); %AA22
                   PorMacroT = zeros(1,1); %AA22 
               end
       end      
       
   end
   
   % DESPLAZAMIENTO IMPUESTO
   vDeltaFixTemp = vfix;
   Du_step_new = zeros(ndoft,1);
   
   % FUERZA EXTERNA IMPUESTA
   if (CONTROL_STRAT == 1)
%       Fext = Fext + f*psi_value; %AA
      Fext = f*psi_value; %AA
      vDeltaFixTemp(doffCondCte) = vfix(doffCondCte)*(psi_value - psi_value_old);
      u(doff) = u(doff) + m_InvCRR*vDeltaFixTemp(doff);
      Du_step_new(doff) = m_InvCRR*vDeltaFixTemp(doff);
   else
      Fext = Fext + f; %AA
      e_VG.vfix =  vfix;
      e_VG.vfix_doff =  m_InvCRR*vfix(doff);
      Du_step_new(doff) =  e_VG.vfix_doff*e_VG.lambda - u(doff) ;
      u(doff) =  e_VG.vfix_doff*e_VG.lambda;
   end
   
   ticIDNewt = tic;
   
   % ESQUEMA DE NEWTON-RAPHSON
   if e_VG.protype==0 % Caso solido
       [u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT,lambda,o_Par] = newton_raphson(xx,m_LinCond,...
           dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,DefMacro,e_VG);
   elseif e_VG.protype==1 % Caso bifasico
       %##########################################################################################
       [u,c_GdlCond,Fint2,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT,lambda,o_Par] = newton_raphsonPM(xx,m_LinCond,...
           dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint2,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
           DefMacro,GradPorMacro,GradPorMacro_up_n,PorMacro,DefMacro_new,GradPorMacro_up_new,PorMacro_new,e_VG);
       %##########################################################################################
   end
   
   %##########################################################################################
   disp('=================================================================');
   fprintf('Tiempo del Newton: %f\n',toc(ticIDNewt));
   disp('=================================================================');
   
   c_ParTT{istep} = o_Par;
   
   % IMPRESION DE RESULTADOS
   index_print = rem(istep,e_VG.IRES);
   if (index_print == 0)
      if mod(istep,postpro_impre_step)==0
          % DESPLAZAMIENTO TOTAL DE LA MICRO-CELDA
          if e_VG.protype==0 %AA
              uTotal = [DefMacroT(1),DefMacroT(4)/2;DefMacroT(4)/2,DefMacroT(2)]*xx(:,1:2)'...
                  +reshape(u,ndn,[]);
          elseif e_VG.protype==1 %AA
              ud=u(e_VG.pos_dG);
              ndn_d=e_VG.ndn_d;
              udTotal = [DefMacroT(1),DefMacroT(4)/2;DefMacroT(4)/2,DefMacroT(2)]*xx(:,1:2)'...
                  +reshape(ud,ndn_d,[]); 
              
            %###################################################################################
             u = f_porp_nint(u,e_DatSet,e_VG); %AA: obtiene los valores de poropresiones en los nodos internos
            %###################################################################################
             up=u(e_VG.pos_pG);
             ndn_p=e_VG.ndn_p;
             upTotal = PorMacroT+GradPorMacroT'*xx(:,1:2)'+reshape(up,ndn_p,[]);
             
         end %AA
         
         %BORRAR
% % % % %          if e_VG.protype==1 %AA
% % % % %              %###################################################################################
% % % % % %              u = f_porp_nint(u,e_DatSet,e_VG); %AA: obtiene los valores de poropresiones en los nodos internos
% % % % %             %###################################################################################
% % % % %              up=u(e_VG.pos_pG);
% % % % %              ndn_p=e_VG.ndn_p;
% % % % %              upTotal = PorMacroT+GradPorMacroT'*xx(:,1:2)'+reshape(up,ndn_p,[]);
% % % % %          end %protype
         
         e_VG.epsilon_Macro = epsilon_Macro;  % AA: Agrego en la estrutura de e_VG porque ingresa en matlab2gid
         e_VG.sigmaHomog = sigmaHomog;  % AA: Agrego en la estrutura de e_VG porque ingresa en matlab2gid
         %AA solo me son utiles para cuando ingresa en f_SaveDatGraf a traves de f_PostProcME
         %$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#
         %$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#
         %$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#
% % % %          e_VG.Dtime = e_VG.theta*time;
         %$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#
         %$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#
         %$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#
         if e_VG.protype==0 %AA
             matlab2gid_res(istep,in,u,c_GdlCond,e_DatSet,e_VarEst_new,e_VarAux,DefMacro,uTotal,...
                 ELOCCalc,e_VG)
         elseif e_VG.protype==1 %AA
             Fext2 = Fext2 + Fext;
             matlab2gid_res_Bif(istep,in,u,c_GdlCond,e_DatSet,e_VarEst_new,e_VarAux,...
                 DefMacro,GradPorMacro,udTotal,upTotal,Fext2,Fint2,ELOCCalc,e_VG)
         end
         
     %###################################################################################
         if e_VG.protype==1 %AA
             u(3*e_VG.in_int) = 0.0; %AA: anulo las poropresiones en los lados (nodos internos)
         end %protype  
    %###################################################################################
      end
   end
   
   if e_VG.isMICRO.MICRO && ~isempty(e_VG.Snap)
       Snapshots = SnapshotSave(Snapshots,istep,e_VarAux,e_VarEst_new,e_DatSet,e_VG,nglT) ;
    end
   %Operaciones Constitutivas despues de la convergencia del Newton
   %Se coloca al final de todo, despues de la impresion de los resultados para que se imprima los
   %valores con los datos que se utilizaron para obtenerlos. Por ejemplo, que las tensiones que se
   %grafica se corresponde con la normal indicada, y no la que se podria obtener del analisis dentro
   %de f_OperConst.
   %#############################################################################################
   [e_DatSet,e_VarEst_new,e_VarAux,sigmaHomog,e_VG] = ...
      f_OperPosConv(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,c_CT,KT,m_LinCond,dofl,doff,e_VG);
   %#############################################################################################
  %AA: FALTA REVISAR PARA QUE APAREZCA eltype=16 y conshyp=14
   
   % Almacenamiento de datos para graficos X-Y
   f_SaveDatGraf(u,c_GdlCond,Fint2,e_VarEst_new,e_VarAux,e_DatSet,m_SetElem,sigmaHomog,epsilon_Macro,e_VG) ;%AA: funciona correctamente 
   
   % ACTUALIZACION DE VARIABLES
   psi_value_old      = psi_value;
   Du_step_old        = Du_step_new;
   e_VarEst_old       = e_VarEst_new;
   if (CONTROL_STRAT == 4)
      e_VG.lambda        = lambda;
   end
   
   if ~mod(istep,deltaPasoSave)
       istepSave = istep;
       save(fullfile(dirFile,'PasoSalvado'))
       %Para que soporte archivos de mas de 1 Gb.
       %save(fullfile(dirFile,'PasoSalvado'),'-v7.3')
   end
   
   fprintf('FIN DE PASO: %-3d (tiempo: %f)\n',istep,toc(ticStep));
   disp('*******************************************************************')
end

%% Almacenamiento de las matrices de Snapshots

if e_VG.isMICRO.MICRO && ~isempty(e_VG.Snap)   
    save([dir '/SNAPSHOTS_' nomb '.mat']  ,'-struct', 'Snapshots');
end







