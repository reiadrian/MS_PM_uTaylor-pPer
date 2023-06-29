function [kuu,kbu_SDA,kbu_MSL,kbbIN_SDA,kbbIN_MSL,fint,res_beta_SDA,...
   res_beta_MSL,sigma_new,hvar_new,eps_new,m_TensorTang,CeBif,...
   ind_state_new,ind_ActState_new,ksd,vectVHElem_new,sigmaTilde_new,m_d_TraccNew] =...
   f_MatElem_quad_q23_SDA_STAFE...
   (delta_u,eps_old,Dbeta_e_SDA,Dbeta_e_MSL,...
   m_VarAuxPG,BI_n,leq_elem,Phi_Grad,n_tensor,CeBif,sigma_old,...
   hvar_old,e_DatElemSet,e_DatMatSet,m_Be,m_DetJe,ind_state_old,ind_ActState_old,...
   e_VG,vectVHElem_old,CPI_n,sigmaTilde_old,volElem)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% factores
%  BI_n           cond. bifurcacion ::: 0 no bifurco; 1 bifurco
%  SL_i           loading condition ::: > 0 carga ; =0 no carga
%  ind_ActState   indice del estado  :::   : 0, 1 o 2
%  CPI            num. lados cortados por el crack path ::: 1,2,3,4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TWO_SCALE      = 0 ;  % comentado para la unificacion de codigos aeh ago_2013

ntens = e_VG.ntens;
ndime = e_VG.ndime;
%
dofpe = e_DatElemSet.dofpe;
npg = e_DatElemSet.npg;
wg = e_DatElemSet.wg;

% Propiedades materiales
sihvarpg = e_DatMatSet.sihvarpg;
%siavarpg = e_DatMatSet.siavarpg;
conshyp  = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;


% Inicializaciones
% ****************
kuu = zeros(dofpe,dofpe);
fint = zeros(dofpe,1);

% SDA ELEMENT
kbu_SDA      = zeros(ndime,dofpe);
kbbIN_SDA    = zeros(ndime,ndime);
res_beta_SDA = zeros(ndime,1);

% MSL ELEMENT
kbu_MSL      = zeros(ntens,dofpe);
kbbIN_MSL    = zeros(ntens,ntens);
res_beta_MSL = zeros(ntens,1);

sigma_new     = zeros(ntens,npg);                     % Included Gauss Points (5 & 6) in the center of the finite element
sigmaTilde_new = zeros(ntens,npg);
eps_new       = zeros(ntens,npg);                       % Included Gauss Points (5 & 6) in the center of the finite element
%Btotal        = zeros(ntens,dofpe);

%Para la homogenizacion del tensor tangente
if esImplex
   m_TensorTang = zeros(ntens,ntens,2*npg);
else
   m_TensorTang = zeros(ntens,ntens,npg);
end

% Redimensionado de matrices
sigma_old        = reshape(sigma_old,ntens,[]);
sigmaTilde_old   = reshape(sigmaTilde_old,ntens,[]);
eps_old          = reshape(eps_old,ntens,[]);
hvar_old         = reshape(hvar_old,sihvarpg,[]);
%aux_var          = reshape(m_VarAuxPG,siavarpg,npg);
hvar_new         = f_DefinicionhVar(conshyp,sihvarpg,npg);
vectVHElem_new   = vectVHElem_old;

%ind_state_new    = 0;
%ind_ActState_new = 0;

if conshyp  == 54
   TWO_SCALE = 1;
end


% SDA PARAMETERS: definicion del ancho de banda
% *********************************************
%Cuando todav�a no fue calculado (previo a la bifurcaci�n) leq_elem en el elemento finito, ya que depende
%depende de la normal de bifurcaci�n.
%No se cambia para el caso multiescala ya que no se hace ninguna regularizaci�n previo a la bifurcaci�n en la
%macroescala (no ser�a necesario al estar regularizados los materiales en la microescala).
if leq_elem==0&&~TWO_SCALE
   %area = m_DetJe(1:4)'*wg(1:4);   %Esto est� mal si no tiene espesor unitario
   %    area = volElem/e_DatElemSet.thickness;
   %    h_eq_zero = sqrt(area);
   %    h_eq = h_eq_zero;
   %Este caso sirve principalmente para el instante de primer da�o, donde el elemento en el paso actual da�a
   %pero reci�n al final de este paso detecta bifurcaci�n (donde se conoce las normales de la discontinuidad y
   %se puede calcular el leq_elem). Este espesor de regularizaci�n se utiliza en ese primer paso en el modelo
   %de da�o. Tambi�n servir�a en los casos (raros) monoescala donde da�a el modelo constitutivo pero no se
   %produce bifurcaci�n.
   %Se toma el valor de la matriz m_hReg.
   switch conshyp
      case {10,11,12,13}  %Modelos de da�o
         h_eq = e_DatMatSet.m_hReg(e_VG.iElemSet);
      otherwise
         error('Cuadr�ngulo SDA bilineal: Longitud de regularizaci�n prefijada: Modelo Constitutivo no definido.')
   end
else
   h_eq = leq_elem;
end

if TWO_SCALE
   %Cuando se usa un modelo constitutivo multiescala kinf es el espesor (medio) del dominio de localizaci�n en
   %la microescala. Como no se conoce previo al c�lculo, en el archivo de datos se ingreso este valor en forma
   %aproximada.
   %En la multiescala se utiliza �nicamente para unos c�lculos
   kinf = e_DatElemSet.kinf;
else
   %En el caso de ser monoescala, kinf es un factor que multiplicado al ancho medio del elemento.
   kinf = e_DatElemSet.kinf*h_eq;
end

% Parametros iniciales para el elemento finito.
switch ind_ActState_old
   case 0
      gamma = 0;
      tau = 1;
      ksd = h_eq;
   case 1
      gamma = 1;
      tau = 0;
      ksd = h_eq;
   case 2
      gamma = 1;
      tau = 0;
      ksd = kinf;
   otherwise
      error('Cuadr�ngulo SDA bilineal: Estado actual del elemento no definido.')
end



if TWO_SCALE == 1
   %    if ksd>10    % VERIFICAR PORQUE Leq SI PUEDE SER >10 DEPENDIENDO DEL GROSOR DEL MALLADO
   %       ksd=kinf;
   %    end
   facNormMicro = hvar_old(6).facNormMicro;
   %facNormMicro = 1;
else
   Eprop_qSD.ksd=ksd;
   Eprop_qSD.gamma=gamma;
   %En el problema de multiescala se realiz� esta correci�n en duro, as� que para que siga funcionando el
   %monoescala se impone que sea 1.
   facNormMicro = 1;
end

m_pesoPG = m_DetJe.*wg;

% B MATRIX ELEMENTAL INTEGRATION
% ******************************
if gamma == 1 % PESOS DE GAUSS INTEGRACION MATRIZ K
   m_pesoPG(5) = volElem;
   m_pesoPG(6) = 0;
end

% *************************************************************
% *************************************************************
% *************************************************************
% *************************************************************
% REGULAR & SINGULAR GAUSS POINT (GP 6)
% STRAIN & STRESS TENSORS ON THE 5th GAUSS POINT (UNLOAD POINT)
% *************************************************************
% *************************************************************
B = m_Be(:,:,6);
m_DefMACRO = B*delta_u;
if ind_ActState_old==2
   % SDA ELEMENTeeee
   %En algunos casos gamma puede ser un valor intermedio.
   m_DefMACRO = m_DefMACRO+gamma*(-Phi_Grad)*Dbeta_e_SDA/facNormMicro;
   %m_DefMACRO = m_DefMACRO+gamma*(-Phi_Grad)*Dbeta_e_SDA;
end

switch conshyp % Modelo constitutivo
   case 11
      if ind_ActState_old==2
         %En algunos casos gamma puede ser un valor intermedio, por eso puede ser que no se quite de esta
         %operaci�n.
         m_DefMACRO = m_DefMACRO + gamma/ksd*n_tensor*Dbeta_e_SDA;
      end
      d_eps_new(:,6) = m_DefMACRO;
      eps_new(:,6) = eps_old(:,6)+d_eps_new(:,6);
      % ONE SCALE MODEL
      e_VG.iPG = 6;
      [ctR,sigma_new(:,6),sigma_new_impl(:,6),hvar_new(:,6)] = ...
         rmap_damage_IMPLEX (...
         eps_new(:,6),hvar_old(:,6),e_DatMatSet,Eprop_qSD,e_VG);
      %
      %% VARIABLE DE INYECCION
      %DISIPACION DE ENERGIA
      %       %Energ�a libre actual
      %       vectVHElem_new(7) = 0.5*((1-hvar_new(1,6))*eps_new(:,6)'*e_DatMatSet.ce*eps_new(:,6));
      %       %Tensi�n media entre dos pasos
      %       sigma_med = 0.5*(sigma_new(:,6)+sigma_old(:,6));
      %       %Potencia interna aplicada entre dos pasos
      %       val_1 = sigma_med'*(eps_new(:,6)-eps_old(:,6));
      %       %Incremento de la energ�a libre entre dos pasos
      %       val_2 =(vectVHElem_new(7)-vectVHElem_old(7));
      %       %Si la potencia interna aplicada es mayor a la energ�a libre recuperable, significa que se disip�
      %       %energ�a.
      %       if val_1>=val_2
      %          %Energ�a disipada por unidad de superficie perpendicular al plano x-y (superficie de fisura)
      %          %hasta el instante actual. Es para comparar con la energ�a de fractura.
      %          %Hay error cuando se produce cambios de ksd (ancho de regularizaci�n) entre pasos??.
      %          vectVHElem_new(1) = vectVHElem_old(1)+ksd*(val_1-val_2);
      %       else
      %          vectVHElem_new(1) = vectVHElem_old(1);
      %       end
      %fprintf('Energ�a total disipada por unidad de superficie %f del elemento %d.\n',...
      %   vectVHElem_new(1),e_VG.iElemNum)
      %
      %�rea bajo la curva q-r del modelo de da�o multiplicado por el espesor de regularizaci�n.
      %Cuando no evoluciona el da�o, es cero el Delta-r impl�cito (hvar_n1(3,6)) , donde el hvar_n1(6,6) es el
      %q actual impl�cito.
      vectVHElem_new(1) = vectVHElem_old(1)+ksd*hvar_new(6,6)*hvar_new(3,6);
      %Se usa siempre kinf para el c�lculo de la energ�a disipada.
      %vectVHElem_new(1) = vectVHElem_old(1)+kinf*hvar_new(6,6)*hvar_new(3,6);
      %       fprintf('�rea bajo la curva q-r por el espesor de regularizaci�n %f del elemento %d.\n',...
      %          vectVHElem_new(1),e_VG.iElemNum)
      %% Variable de suavizado del crack path field.
      vectVHElem_new(4) = hvar_new(4,6);
      %% Criterio de descarga del elemento
      %Se utiliza el incremento impl�cito de r.
      %vectVHElem_new(5) = hvar_new(3,6);
      %Se utiliza el incremento expl�cito de r.
      vectVHElem_new(5) = hvar_new(8,6);
      %% Variable r de trial
      vectVHElem_new(9) = hvar_new(9,6);
      % *******************************************
      delta_sigma(:,6) = sigma_new(:,6)-sigma_old(:,6);
      %Se calcula las tracciones totales.
      m_d_TraccNew = n_tensor'*sigma_new(:,6);
      %
      if ind_ActState_old==2
         %if ind_state_old==2
         if esImplex
            m_CTHomog_ub = n_tensor'*ctR.Implex;
            m_CTHomog_bb = n_tensor'*ctR.Implex*n_tensor/ksd;
            ctR.Implex(1:2,1:2)   = m_CTHomog_bb;
            ctR.Implex(3:4,:) = m_CTHomog_ub;
         else
            m_CTHomog_ub = n_tensor'*ctR;
            m_CTHomog_bb = n_tensor'*ctR*n_tensor/ksd;
            ctR.Impli(1:2,1:2)   = m_CTHomog_bb;
            ctR.Impli(3:4,:)     = m_CTHomog_ub;
         end
      end
   case 54   % TWO SCALE MODEL  SANTA FE PG 6
      e_VG.iPG   = 6;
      e_VG.elast = 0;
      d_eps_new(:,6) = m_DefMACRO;
      %
      m_IDefR = m_DefMACRO;
      m_IDefS = Dbeta_e_SDA;
      %
      %La longitud de regularizaci�n kinf ingresada se utiliza �nicamente en la funci�n
      %Regularization_parameters_Micro para el c�lculo de la energ�a disipada y la inyecci�n de la SD, que se
      %calcula solo para el PG 6. Como se quiere tener un par�metro de referencia de cu�nto energ�a de
      %fractura se consumi� se debe ingresar el espesor de la banda en la microescala, que tiene que tener
      %relaci�n con kinf ingresado en la macroescala (una media de los espesores micro en la banda en la
      %microescala).
      %El ksd se ingresa para obtener, en el estado 1, un resultado aproximadamente objetivo con el tama�o de
      %la microcelda.
      %En m_d_TraccNew se lleva las tracciones totales ya que ahora se hace equilibrio de las fuerzas internas
      %totales.
      [ctR,sigma_new(:,6),hvar_new(:,6),vectVHElem_new,m_d_TraccNew] =...
         f_RMap_MEStafe(m_IDefR,m_IDefS,hvar_old(:,6),...
         e_DatMatSet,BI_n,ind_ActState_old,...
         e_VG,vectVHElem_old,ksd,kinf);
      % *******************************************
      delta_sigma(:,6) = sigma_new(:,6)-sigma_old(:,6);
end

if esImplex
   c_tilde6              = ctR.Implex;
   m_TensorTang(:,:,6)   = c_tilde6 ;
   m_TensorTang(:,:,12)  = ctR.Impli;
else
   c_tilde6              = ctR;
   m_TensorTang(:,:,6)   = c_tilde6 ;
end

% if ind_state_old<=0
%    if esImplex
%       CeBif     = m_TensorTang(:,:,6);
%    else
%       if e_VG.istep==1
%          CeBif     = m_TensorTang(:,:,6);
%       end
%    end
% end
% REVISION DE LOS PARAMETROS DE INYECCION DE LA SD
if (BI_n == 0) % SD Range
   %    if ~isnan(e_DatElemSet.gfvRef)
   %       if TWO_SCALE
   %          %Seleccionar en duro el set material en la microescala donde se extrae la energ�a de fractura de
   %          %referencia.
   %          iSet = 3;
   %          FRAC_ENER = e_DatMatSet.e_DatSet(iSet).e_DatMat.gfv;
   %       else
   %          FRAC_ENER = e_DatMatSet.gfv;
   %       end
   %    else
   FRAC_ENER = e_DatElemSet.gfvRef;
   %    end
   facIny = e_DatElemSet.facIny;
   vectVHElem_new(2) = facIny*FRAC_ENER;  % FRACTURE ENERGY CRITERIA
   vectVHElem_new(3) =  100000.0*facIny*vectVHElem_new(2); % RANGO NO USADO (POR SI SE QUIERE USAR 3 TIPOS DE ELEMENTOS FINITOS)
else
   vectVHElem_new(2) = vectVHElem_old(2)   ; %injection  % en el PG=6 (singular) almaceno las variables elementales
   vectVHElem_new(3) = vectVHElem_old(3)   ; %Posibilidad de un tercer tipo de EF
end
% *************************************************************
% *************************************************************
% *************************************************************
% *************************************************************
% REGULAR GAUSS POINT (GP 5)
% STRAIN & STRESS TENSORS ON THE 5th GAUSS POINT (UNLOAD POINT)
% *************************************************************
% *************************************************************
%%%% if ind_ActState_new == 0 && BI_n == 0
%S�lo se copia hasta que se activa alguna vez la SD (se lleg� por primera vez al estado 2), luego se asume que
%evolucionan por separado.
%%
%if ind_ActState_old<=1
if ind_state_old<=1
   
   %No ser�a necesario hacer esta copia, excepto para imprimir, durante la etapa que usa solo 4 puntos de
   %gauss (ind_ActState_old==0) y durante etapa de 1 solo puntos de gauss hasta la activaci�n de la SD
   %(ind_ActState_old==1).
   hvar_new(:,5)   = hvar_new(:,6);
   sigma_new(:,5)  = sigma_new(:,6);
   %Como en el caso del elemento 54 no se calcula una deformaci�n singular en el PG6 ya que no se utiliza, y 
   %deja de tener sentido en el contexto multiescala despu�s que se inyecta, se calcula para el PG5.
   switch conshyp % Modelo constitutivo
      case 11
         eps_new(:,5) = eps_new(:,6);
      case 54
         %Hasta la inyecci�n (que comienza evolucionar Beta) m_DefMACRO se corresponde con la que se calcula
         %para el PG6.
         d_eps_new(:,5) = d_eps_new(:,6);
         eps_new(:,5) = eps_old(:,5)+d_eps_new(:,5);
      otherwise
         error('Punto de gauss 5: Actualizaci�n de variables: Elemento Finito no definido.')   
   end      
   %      sigma_new_impl(:,5) = sigma_new_impl(:,6);
   if esImplex
      c_tilde5            = c_tilde6;
      m_TensorTang(:,:,5) = c_tilde5 ;
      m_TensorTang(:,:,11)= m_TensorTang(:,:,12);
   else
      c_tilde5            = c_tilde6;
      m_TensorTang(:,:,5) = c_tilde5;
   end
   delta_sigma(:,5) = delta_sigma(:,6);
   
   if (gamma == 1)
      kuu = kuu  + (B'* c_tilde5 *B)*m_pesoPG(5);
   end
   
else
   
   B = m_Be(:,:,5);
   m_DefMACRO = B*delta_u;
   %m_DefMACRO = m_DefMACRO-gamma*Phi_Grad*Dbeta_e_SDA;
   m_DefMACRO = m_DefMACRO-gamma*Phi_Grad*Dbeta_e_SDA/facNormMicro;
   
   d_eps_new(:,5) = m_DefMACRO;
   eps_new(:,5) = eps_old(:,5) + d_eps_new(:,5) ;
   switch conshyp % Modelo constitutivo
      case 11   % modelo da�o monoescala %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         e_VG.iPG = 5;
         if ind_ActState_old<=1
            [ctR,sigma_new(:,5),sigma_new_impl(:,5),hvar_new(:,5)] = ...
               rmap_damage_IMPLEX (eps_new(:,5),hvar_old(:,5),e_DatMatSet,Eprop_qSD,e_VG);
         else
            [ctR,sigma_new(:,5),sigma_new_impl(:,5),hvar_new(:,5)]  = ...
               rmap_damage_elastic(eps_new(:,5),hvar_old(:,5),e_DatMatSet);
         end
      case 54   %Modelo multiescala cohesivo SANTA FE PG 5
         e_VG.iPG = 5;
         %%%%%%%%%%%%%%%%%% VIEJO    %%%%%%%%%%%%%%%%%%%
%                   if ind_ActState_old<=1
%                      e_VG.elast = 0;
%                   else
%                      e_VG.elast = 1;
%                   end
%                   condBifR = 0;
%                   m_IDefR = d_eps_new(:,5);
%                   m_IDefS = zeros(ntens,1);
%                   %
%                   [ctR,sigma_new(:,5),hvar_new(:,5)] =...
%                      f_RMap_MEStafe(m_IDefR,m_IDefS,hvar_old(:,5),...
%                      e_DatMatSet,condBifR,ind_ActState_old,...
%                      e_VG,vectVHElem_old,ksd);
         %
         %%%%%%%%%%%%%%%%%% NUEVO   %%%%%%%%%%%%%%%%
         if ind_ActState_old<=1
            e_VG.elast = 0;
            condBifR = 0;
            m_IDefR = d_eps_new(:,5);
            m_IDefS = zeros(ntens,1);
            %
            [ctR,sigma_new(:,5),hvar_new(:,5)] =...
               f_RMap_MEStafe(m_IDefR,m_IDefS,hvar_old(:,5),...
               e_DatMatSet,condBifR,ind_ActState_old,...
               e_VG,vectVHElem_old,ksd);
         else
            sigma_new(:,5) = CeBif*d_eps_new(:,5)+sigma_old(:,5);
            hvar_new(:,5) = hvar_old(:,5);
            if esImplex
               ctR.Implex = CeBif;
               ctR.Impli = CeBif;
            else
               ctR = CeBif;
            end
         end
         %%%%%%%%%%%%%%%%%% NUEVO   %%%%%%%%%%%%%%%%
   end
   
   if esImplex
      c_tilde5              = ctR.Implex;
      m_TensorTang(:,:,5)   = c_tilde5 ;
      m_TensorTang(:,:,11)  = ctR.Impli;
   else
      c_tilde5              = ctR;
      m_TensorTang(:,:,5)   = c_tilde5 ;
   end
   
   delta_sigma(:,5) = sigma_new(:,5)-sigma_old(:,5);
   if gamma==1 % SD Range
      kuu = kuu+(B'*c_tilde5*B)*m_pesoPG(5);
   end
   
end
% *********************************************************************
% *********************************************************************
% **************  Puntos de gauss regulares        ********************
% *********************************************************************

% STRAIN & TENSION TENSORS ON THE STANDARD (REGULAR) GAUSS POINTS (1-4)
% *********************************************************************
for iPG = 1:(npg-2) % swap the initial four gauss points - regular points
   B = m_Be(:,:,iPG);
   m_DefMACRO  = B*delta_u;
   if ind_ActState_old==2
      % SDA_ELEMENT
      %m_DefMACRO = m_DefMACRO - gamma*Phi_Grad*Dbeta_e_SDA;
      m_DefMACRO = m_DefMACRO - gamma*Phi_Grad*Dbeta_e_SDA/facNormMicro;
   end
   
   d_eps_new(:,iPG)  = m_DefMACRO;
   eps_new(:,iPG)    = eps_old(:,iPG) + d_eps_new(:,iPG);
   switch conshyp % Modelo constitutivo
      %case 12
      case 11
         %------------------
         %Los PG de 1 a 4 siempre se regularizan con el tama�o del elemento. Esto solo influye en el caso de
         %haber descarga y posterior recarga.
         Eprop_qSD_4PG.ksd = h_eq;
         %------------------
         e_VG.iPG = iPG;
         [ctR,sigma_new(:,iPG),sigma_new_impl(:,iPG),hvar_new(:,iPG)] = ...
            rmap_damage_IMPLEX (...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,Eprop_qSD_4PG,e_VG);
      case 54   % TWO SCALE MODEL  PG : iPG
         e_VG.iPG   = iPG;
         e_VG.elast = 0;
         condBifR   = 0;
         m_IDefR    = d_eps_new(:,iPG);
         m_IDefS    = zeros(ntens,1);
         %
         %Los PG de 1 a 4 siempre se regularizan con el tama�o del elemento. Esto solo influye en el caso de
         %haber descarga y posterior recarga.
         [ctR,sigma_new(:,iPG),hvar_new(:,iPG),~, ~] =...
            f_RMap_MEStafe(m_IDefR,m_IDefS,hvar_old(:,iPG),...
            e_DatMatSet,condBifR,ind_ActState_old,...
            e_VG,vectVHElem_old,h_eq);
   end
   if esImplex
      c_tilde                   = ctR.Implex;
      m_TensorTang(:,:,iPG)     = c_tilde ;
      m_TensorTang(:,:,iPG+npg) = ctR.Impli;
   else
      c_tilde                   = ctR;
      m_TensorTang(:,:,iPG)     = c_tilde ;
   end
   
   if (gamma == 0) % SD Range
      kuu = kuu  + tau*(B'* c_tilde *B)*m_pesoPG(iPG);
   end
   
   %Sebastian: Considerando siempre es psi nulo, como est� actualmente.
   %       delta_sigma(:,iPG) = (1-gamma)*(tau*(sigma_new(:,iPG)-sigma_old(:,iPG))+...
   %           (1-tau)*(sigma_new(:,6)-sigma_old(:,6)))...
   %           + ...
   %           gamma*(psi*(sigma_new(:,6)-sigma_old(:,6)) + (1-psi)*(sigma_new(:,5)-sigma_old(:,5)));
   delta_sigma(:,iPG) = (1-gamma)*(tau*(sigma_new(:,iPG)-sigma_old(:,iPG))+...
      (1-tau)*(sigma_new(:,6)-sigma_old(:,6)))...
      + gamma*(sigma_new(:,5)-sigma_old(:,5));
   
   sigmaTilde_new(:,iPG) = sigmaTilde_old(:,iPG) + delta_sigma(:,iPG);
   
   % Calculo de fint
   % ***************
   % fint = fint+B'*delta_sigma(:,iPG)*m_pesoPG(iPG);
   fint = fint+B'*sigmaTilde_new(:,iPG)*m_pesoPG(iPG);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if gamma==1   %Sebastian: No es necesario la verificaci�n con gamma==1 si despu�s se verifica s�lo con
%ind_ActState_old (distinto si se verificara con ind_ActState_old==1, que tambi�n tiene gamma==1)
if ind_ActState_old==2
   %La longitud s solo es un factor de estabilizaci�n, para que las matrices ensambladas tengan mejor
   %n�mero de condici�n (el resultado es independiente de este valor).
   s = volElem/h_eq;
   % SDA ELEMENT
   % ***********
   kub_SDA = -m_Be(:,:,5)'*c_tilde5*Phi_Grad/facNormMicro*m_pesoPG(5) ;
   m_Dctn  = c_tilde6(ndime+1:2*ndime,:)+n_tensor'*(-c_tilde5);
   kbu_SDA = s*m_Dctn*m_Be(:,:,5) ;
   %
   rCondMat = rcond(s*(c_tilde6(1:ndime, 1:ndime)-m_Dctn*Phi_Grad)/facNormMicro);
   if rCondMat<eps(1)*100||isnan(rCondMat)
      fprintf('Matriz de discontinuidad, del elemento %d, mal condicionada (rcond=%g).\nMatriz:\n',...
         e_VG.iElemNum,rCondMat)
   end
   %
   kbbIN_SDA = (s*(c_tilde6(1:ndime, 1:ndime)-m_Dctn*Phi_Grad)/facNormMicro)\eye(ndime);
   
   % TRACTION CONTINUITY OVER THE LOCALIZATION BAND
   % **********************************************
   res_beta_SDA = s*(m_d_TraccNew-n_tensor'*sigma_new(:,5)) ;
   % K MATRIX AND RESIDUAL CONDENSED FORCE
   % *************************************
   m_k1 = kub_SDA*kbbIN_SDA;
   kuu = kuu-m_k1*kbu_SDA;
   
   fint = fint - m_k1*res_beta_SDA;
end
%end

%Actualizaci�n de variables de estados del elemento (en que etapa se etapa se encuentra el SD). Se asume una
%actualizaci�n expl�cita del estado.
%Todas estos cambios de estados se cambia en la funci�n posconvergencia.
% if BI_n==1
%
%    switch conshyp
%       case {4,11,44,53,54}
%          ro = vectVHElem_new(1) ;       %  hvar_old(15,6);
%          ro_inj = vectVHElem_new(2) ;   %  hvar_old(8,6);
%       otherwise
%          error('Cuadr�ngulo SDA bilineal: Factores ro: Modelo Constitutivo no definido.')
%    end
%    SLI_n = vectVHElem_new(5);
%
%    % FROM (STATE==0) TO (STATE==1 || STATE==2)
%    if ind_ActState_old==0
%
%       %No existe la posibilidad que sea CPI_n==3???? En un cuadr�ngulo parece imposible, como el caso de
%       %CPI_n==1.
%       % CASE 1
%       %if ((RLI_n>0) && (ro<ro_inj)) || (((RLI_n>0) && (ro>ro_inj)) && (CPI_n==0 || CPI_n==4))
%       if ((SLI_n>0) && (ro<ro_inj)) || (((SLI_n>0) && (ro>ro_inj)) && (CPI_n==0 || CPI_n==4))
%          ind_state_new = 1;
%          ind_ActState_new = 1;
%       % CASE 2
%       %elseif (((SLI_n>0) && (ro>ro_inj)) || ((RLI_n>0) && (ro>ro_inj))) && (CPI_n==2)
%       elseif (((SLI_n>0) && (ro>ro_inj)) || ((SLI_n>0) && (ro>ro_inj))) && (CPI_n==2)
%          ind_state_new = 2;
%          ind_ActState_new = 2;
%       else
%          ind_ActState_new = ind_ActState_old;
%          ind_state_new = ind_state_old;
%       end
%
%    % FROM (STATE==1) TO (STATE==0 || STATE==2)
%    elseif ind_ActState_old==1
%
%       % CASE 3
%       if (SLI_n<=0)
%          ind_state_new = 1;
%          ind_ActState_new = 0 ;
%       % CASE 4
%       elseif (SLI_n>0) && (ro>ro_inj) && (CPI_n==2)
%          ind_state_new = 2;
%          ind_ActState_new = 2;
%       else
%          ind_ActState_new = ind_ActState_old;
%          ind_state_new = ind_state_old;
%       end
%
%    % FROM (STATE==2) TO (STATE==0 || STATE==1)
%    elseif ind_ActState_old == 2
%
%       % CASE 5
% %       if (SLI_n<=0)
% %          ind_state_new = 2;
% %          ind_ActState_new = 0 ;
% %       % CASE 6
% %       elseif (SLI_n>0) && (CPI_n==0 || CPI_n==4)
% %          ind_state_new = 2;
% %          ind_ActState_new = 1;
% %       else
% %          ind_ActState_new = ind_ActState_old;
% %          ind_state_new = ind_state_old;
% %       end
%       %Si no se quiere que descargue del SD.
%       ind_ActState_new = ind_ActState_old;
%       ind_state_new = ind_state_old;
%
%    end
%
% else
%
ind_ActState_new = ind_ActState_old;
ind_state_new = ind_state_old;
%
% end

%Se ordena las matrices como vectores columnas (es mas rapido usar Mat(:), que reshape(Mat,[],1);).
sigma_new = sigma_new(:);
sigmaTilde_new = sigmaTilde_new(:);
eps_new = eps_new(:);
hvar_new = hvar_new(:);
%   sigma_tilde_new = sigma_tilde_new(:);

end