%Determina operadores tangentes homogeneizados, tensiones homogeneizadas, cantidad de masa
%del fuido homogeneizada, velocidades homogeneizadas y otras variables
function [e_TanOp,sigmaE_new,sigmaT_new,mflu_new,velflu_sta,velflu_total,velflu_new,...
    hvar_newMacro] = f_RMap_ME_Bifase(eps_new,phi_new_dup,phi_new_up_n,...
                        p_new,eps_n1,phi_new_up_n1,p_n1,hvar_oldMacro,...
                        e_DatMatSetMacro,e_VGMacro)

% if e_VGMacro.iElemNum==6&&e_VGMacro.iPG==1
%     a=1;
% elseif e_VGMacro.iElemNum==5&&e_VGMacro.iPG==1
%     a=1;
% elseif e_VGMacro.iElemNum==4&&e_VGMacro.iPG==1
%     a=1;
% elseif e_VGMacro.iElemNum==3&&e_VGMacro.iPG==1
%     a=1;
% elseif e_VGMacro.iElemNum==2&&e_VGMacro.iPG==1
%     a=1;
% elseif e_VGMacro.iElemNum==1&&e_VGMacro.iPG==1
%     a=1;
% end
   
%Se recupera variables micro
xx = e_DatMatSetMacro.xx;
omegaMicro_d = e_DatMatSetMacro.omegaMicro_d;
e_DatSet = e_DatMatSetMacro.e_DatSet;
e_VG = e_DatMatSetMacro.e_VG;
esImplexMacro = e_DatMatSetMacro.esImplex;
e_VarEst_old = hvar_oldMacro.e_VarEst;
% Vector de incognitas al paso de tiempo previo "n"
u = hvar_oldMacro.u;
c_GdlCond = hvar_oldMacro.c_GdlCond;
Fint         = hvar_oldMacro.Fint;
m_LinCond    = hvar_oldMacro.m_LinCond;
doff         = hvar_oldMacro.doff;
dofl         = hvar_oldMacro.dofl;

e_VarAux = hvar_oldMacro.e_VarAux;

% VARIABLES GLOBALES
nElemTot = e_VG.nElem;
ntens = e_VG.ntens;
ndoft = e_VG.ndoft;

% INICIALIZACION DE VARIABLES
% Vector de incrementos de desplazamientos en el paso de tiempo previo
Du_step_old = zeros(ndoft,1); 

% Vector de incrementos de fuerzas externas micro-escala
Fext = zeros(ndoft,1);

%Por si es necesario el paso tiempo a nivel micro.
e_VG.istep = e_VGMacro.istep;
e_VG.Dtime = e_VGMacro.Dtime;
e_VG.Dtime2 = e_VGMacro.Dtime2;
% Guardo conshyp MACRO para calcular velflu en f_MatElem_BifaseMulTSc
e_VG.conshypMacro = e_VGMacro.conshyp;

%Se guarda que se quiere que el modelo constitutivo se comporte elasticamente.
e_VG.elast = e_VGMacro.elast;

%Para impresion y debug se guarda la iteracion macro.
e_VG.iterMacro = e_VGMacro.iter;
%Para imprension se guarda el numero de elemento macro y numero de PG macro.
e_VG.iElemNumMacro = e_VGMacro.iElemNum;
e_VG.iPGMacro = e_VGMacro.iPG;
e_VG.SharedFS = e_VGMacro.SharedFS; 

%Se guarda los datos del matlabpool (se utiliza para imprimir de pasos y iteraciones a nivel micro que no
%convergieron).
e_VG.nLab = e_VGMacro.nLab; 
e_VG.tipoMPool = e_VGMacro.tipoMPool; 

% DESPLAZAMIENTO IMPUESTO
%No seria necesario estas operaciones porque vfix en todos los modelos clasicos de las
%formulaciones multiescala son nulos. Ademas que habria que interpretar como realizar el delta
%psi_value a nivel micro, ya deberia se corresponder con el incremento de tiempo a nivel macro,
%pero lo que implicaria que en cada paso de tiempo se esta resolviendo un RVE distinto.

% VARIABLE PRIMITVA macro aplicada en la Celda unitaria.
%Se aplica en forma uniforme en todo dominio.

%Deformacion macro por elemento y por punto de gauss, dividida en sets.
%Delta de macro-deformacion
c_DefMacro = arrayfun(@(x)repmat(eps_new,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false); 
%Delta de macro gradiente de poro presiones
c_GradPorMacro_dup = arrayfun(@(x)repmat(phi_new_dup,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
%Macro gradiente de poropresiones en el tiempo "n"
c_GradPorMacro_up_n = arrayfun(@(x)repmat(phi_new_up_n,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
%Delta de macro-poro presiones
c_PorMacro = arrayfun(@(x)repmat(p_new,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);


%Macro-deformacion en el tiempo "n+1"
c_DefMacro_new = arrayfun(@(x)repmat(eps_n1,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false); 
%Macro gradiente de poropresiones en el tiempo "n+1"
c_GradPorMacro_up_new = arrayfun(@(x)repmat(phi_new_up_n1,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
%Macro-poro presiones en el tiempo "n+1"
c_PorMacro_new = arrayfun(@(x)repmat(p_n1,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);


%Deformacion macro por elemento, dividad en sets.
%    c_DefMacro = arrayfun(@(x)repmat(eps_new,[1,x.nElem]),e_DatSet,'UniformOutput',false);
%    c_GradPorMacro = arrayfun(@(x)repmat(phi_new,[1,x.nElem]),e_DatSet,'UniformOutput',false);
%    c_PorMacro = arrayfun(@(x)repmat(p_new,[1,x.nElem]),e_DatSet,'UniformOutput',false);


Du_step_new = zeros(ndoft,1);

% ESQUEMA DE NEWTON-RAPHSON
if e_VGMacro.conshyp~=62
    [u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] = newton_raphsonPM(...
    xx,m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,...
    e_VarAux,e_DatSet,c_DefMacro,c_GradPorMacro_dup,c_GradPorMacro_up_n,c_PorMacro,...
    c_DefMacro_new,c_GradPorMacro_up_new,c_PorMacro_new,e_VG);
elseif  e_VGMacro.conshyp==62
    [u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] = newton_raphson_analyticPM(...
    xx,m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,...
    e_VarAux,e_DatSet,c_DefMacro,c_GradPorMacro_dup,c_GradPorMacro_up_n,c_PorMacro,...
    c_DefMacro_new,c_GradPorMacro_up_new,c_PorMacro_new,e_VG,omegaMicro_d);
end


% OPERADORES TANGENTES HOMOGENEIZADOS
e_TanOp = f_OpTangHomBifase(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro_d,...
true(nElemTot,1),true(nElemTot,1),xx,e_VG);

%Se asume que no se realiza analisis de bifurcacion con el tensor tangente constitutivo homogeneizado, por
%lo que en el caso ser implex, se devuelve nulo el tensor implicito homogeneizado.
% VEEEEER!!!!!!!!!!!!!!
if esImplexMacro
    m_CTHomog = struct('Implex',m_CTHomog,'Impli',zeros(ntens,ntens));
end
%##############################################################################################
%CALCULO DE VARIABLES DUALES HOMOGENEIZADAS
%Delta de tension efectiva homogeneizada
sigmaE_new = f_HomogArea({e_VarEst_new.sigmaE},ntens,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
%Delta de tension total homogeneizada
sigmaT_new = f_HomogArea({e_VarEst_new.sigmaT},ntens,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
%Delta de cantidad de masa del fluido homogeneizada
mflu_new = f_HomogArea({e_VarEst_new.mflu},1,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
%Velocidad de filtracion homogeneizada al tiempo "n+theta' 
velflu_sta = f_HomogArea({e_VarEst_new.velflu_sta},2,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
%Velocidad de filtracion al paso de tiempo "n+1" 
velflu_total = f_HomogArea({e_VarEst_new.velflu_total},2,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
%Velocidad de filtracion al paso de tiempo "n+theta" 
velflu_new = f_HomogArea({e_VarEst_new.velflu},2,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
%##############################################################################################

%##############################################################################################
%CALCULO DE MEDIAS VOLUMETRICAS
%Variables que deben tender a cero (dada la teoria)
%Delta del gradiente de micro-poro presiones fluctuantes homogeneizadas
grad_p_HomogFl= f_HomogArea({e_VarEst_new.phi_fluct},2,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
%Delta de micro-deformaciones fluctuantes homogeneizadas
defHomogFl = f_HomogArea({e_VarEst_new.eps_fluct},ntens,omegaMicro_d,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
%Micro-desplazamientos fluctuantes y micro-poro presiones fluctuantes
%homogeneizados al tiempo "n+1"
%Verificacion de la media del desplazamiento y de la poropresion
[m_uMedioFl_d,m_uMedioFl_p] = f_MediaDespCU_Bif(u,omegaMicro_d,e_DatSet,e_VG);

fprintf('Elemento %d: PG %d: ||Delta deformacion fluctuante media||     : %g\n',e_VGMacro.iElemNum,...
e_VGMacro.iPG,norm(defHomogFl))
fprintf('Elemento %d: PG %d: ||Desplazamiento fluctuante medio||    : %g\n',e_VGMacro.iElemNum,...
e_VGMacro.iPG,norm(m_uMedioFl_d))
fprintf('Elemento %d: PG %d: ||Delta Grad(poro presion) fluctuante media||: %g\n',e_VGMacro.iElemNum,...
e_VGMacro.iPG,norm(grad_p_HomogFl)) 
fprintf('Elemento %d: PG %d: ||Poro presion fluctuante media||    : %g\n',e_VGMacro.iElemNum,...
e_VGMacro.iPG,norm(m_uMedioFl_p))

%##############################################################################################
if norm(grad_p_HomogFl)>1*10^(-14) || norm(m_uMedioFl_p)>1*10^(-14) ...
        || norm(defHomogFl)>1*10^(-14) || norm(m_uMedioFl_d)>1*10^(-14)
    error('ALGUNA MEDIA NO SATISFACE SER NULA (NUMERICAMENTE)')
end
%##############################################################################################

%##############################################################################################
%SALIDA DE INFORMACION MICRO-ESCALA
for iSet = 1:e_VG.nSet
    %############################################################################################
    %VARIABLES MICRO-ESCALA EN EL PASO DE TIEMPO "n+1"
    %############################################################################################
    e_VarEst_new(iSet) = f_OutVar_time(e_VarEst_new(iSet),e_VarEst_old(iSet),...
            e_VarEst_new(iSet).eps,e_VarEst_new(iSet).phi,...
            e_VarEst_new(iSet).porpr,e_VarEst_new(iSet).eps_fluct,...
            e_VarEst_new(iSet).phi_fluct,e_VarEst_new(iSet).p_fluct,...
            e_VarEst_new(iSet).sigmaE,e_VarEst_new(iSet).sigmaT,...
            e_VarEst_new(iSet).mflu,e_VarEst_new(iSet).velflu_sta,...
            e_VarEst_new(iSet).velflu_total,e_VarEst_new(iSet).velflu,e_VG,0);
    %############################################################################################
    %############################################################################################
    %VARIABLES MICRO-ESCALA EN EL PASO DE TIEMPO "n+1"
    %############################################################################################
%     e_VarEst_new(iSet) = f_OutVar_time(e_VarEst_new(iSet),e_VarEst_old(iSet),...
%             e_VarEst_new(iSet).eps,e_VarEst_new(iSet).phi,...
%             e_VarEst_new(iSet).porpr,e_VarEst_new(iSet).eps_fluct,...
%             e_VarEst_new(iSet).phi_fluct,e_VarEst_new(iSet).p_fluct,...
%             e_VarEst_new(iSet).sigmaE,e_VarEst_new(iSet).sigmaT,...
%             e_VarEst_new(iSet).mflu,e_VarEst_new(iSet).velflu_sta,...
%             e_VarEst_new(iSet).velflu_dyn,e_VG,1);
    %############################################################################################
    %#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~
    %NO ES PRACTICO. SUPONE NPG IGUAL EN CADA SET Y ME QUEDARIA CON EL ULTIMO
    npg=e_DatMatSetMacro.e_DatSet(iSet).e_DatElem.npg; 
    %#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~#!~
    if e_VGMacro.istep>1
        nElem = e_DatSet(iSet).nElem;
        c_DefMacro(iSet) = mat2cell(cell2mat(c_DefMacro(iSet)) + cell2mat(hvar_oldMacro.c_DefMacro(iSet)),ntens,npg,nElem);
        c_GradPorMacro_dup(iSet) = mat2cell((cell2mat(hvar_oldMacro.c_GradPorMacro_dup(iSet)) + cell2mat(c_GradPorMacro_dup(iSet))),2,npg,nElem);
    end
end
%############################################################################################

%###################################################################################
%Agregue Fext
hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,'c_DefMacro',{c_DefMacro},...
'c_GradPorMacro_dup',{c_GradPorMacro_dup},'c_GradPorMacro_up_n',{c_GradPorMacro_up_n},...
'c_PorMacro',{c_PorMacro},'c_DefMacro_new',{c_DefMacro_new},'c_GradPorMacro_up_new',{c_GradPorMacro_up_new},...
'c_PorMacro_new',{c_PorMacro_new},'Fext',Fext);
 %###################################################################################
end

%%
%Determina los micro-desplazamientos fluctuantes y micro-poro presiones fluctuantes homogeneizados
function [m_uMedio_d,m_uMedio_p] = f_MediaDespCU_Bif(u,omegaMicro,e_DatSet,e_VG)

   nSet = e_VG.nSet;
   %El desplazamiento es un campo nodal, por lo que para intregrar sobre el dominio se lo lleva a los punto de
   %Gauss.
   c_DespElem = cell(nSet,1);
   c_PorpElem = cell(nSet,1);

   %
   for iSet = 1:nSet
      %
      nElem = e_DatSet(iSet).nElem;
      nPG = e_DatSet(iSet).e_DatElem.npg;
      m_DofElem = e_DatSet(iSet).m_DofElem;
      m_FF_d = e_DatSet(iSet).m_FF_d;
      m_FF_p = e_DatSet(iSet).m_FF_p;
      %
      pos_d =  e_DatSet(iSet).e_DatElem.pos_d;
      pos_p =  e_DatSet(iSet).e_DatElem.pos_p;
      m_uElemSet_d = reshape(u(m_DofElem(pos_d,1:nElem)),[],nElem);
      m_uElemSet_p = reshape(u(m_DofElem(pos_p,1:nElem)),[],nElem);
      %
      m_DespPG = zeros(2,nPG,nElem);
      m_PorpPG = zeros(1,nPG,nElem);
      for iElem = 1:nElem
         %squeeze llama reshape, asi que no es mas rapido que esta si se conoce cual es la dimension de la
         %matriz con valor 1.
         %Desplazamientos
         m_DespPG(:,:,iElem) = squeeze(sum(bsxfun(@times,m_FF_d(:,:,:),m_uElemSet_d(:,iElem)'),2));
         %Poropresiones
         m_PorpPG(:,:,iElem) = squeeze(sum(bsxfun(@times,m_FF_p(:,:,:),m_uElemSet_p(:,iElem)'),2));
      end
      %
      c_DespElem{iSet} = m_DespPG;
      c_PorpElem{iSet} = m_PorpPG;
   end
   
   m_uMedio_d = f_HomogArea(c_DespElem,2,omegaMicro,{e_DatSet.m_DetJT_d},e_DatSet,e_VG);
   m_uMedio_p = f_HomogArea(c_PorpElem,1,omegaMicro,{e_DatSet.m_DetJT_p},e_DatSet,e_VG);

end
