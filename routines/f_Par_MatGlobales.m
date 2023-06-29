function [KT,Fint,c_GdlCond,e_VarEst_new,e_VarAux,c_TensorTang,o_Par] = f_Par_MatGlobales(xx,u,du,...
   c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,e_VG)

o_Par = [];
%Para no transferir estas variables que no se utilizan dentro del parfor.
e_VG.smooth_alpha = [];
e_VG.smooth_dalpha = [];
e_VG.smooth_dalpha_old = [];
e_VG.MGlobInv = [];
%[e_DatSet(:).dN_xy] = deal([]);
%[e_DatSet(:).MElemInv] = deal([]);
%[e_DatSet(:).N_vector] = deal([]);

% Recupera variables globales
ntens = e_VG.ntens;
nSet = e_VG.nSet;
ndime = e_VG.ndime;
protype = e_VG.protype;

% Inicializaciones
c_Ke = cell(nSet,1);
c_Fint = cell(nSet,1);
c_Fil = cell(nSet,1);
c_Col = cell(nSet,1);
c_FilFza = cell(nSet,1);
c_TensorTang = cell(nSet,1);

%Loop sobre el set de elementos (tipo y material del elemento)
for iSet = 1:nSet
   
   % Recuperacion de variables
   %(evita llamar desde las estructura de los sets en el bucle de los elementos, ver que solo
   %recupera punteros)
   nElem = e_DatSet(iSet).nElem;
   %conec = e_DatSet(iSet).conec;
   m_DofElem = e_DatSet(iSet).m_DofElem;
   
   %sigma_old = e_VarEst_old(iSet).sigma; AA: meti dentro de if-elseif mas
   %abajo
   eps_old = e_VarEst_old(iSet).eps;
   hvar_old = e_VarEst_old(iSet).hvar;
   m_VarHistElemOld = e_VarEst_old(iSet).VarHistElem;
   aux_var = e_VarAux(iSet).VarAuxGP;
   m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
   %Las funciones de forma en coordenadas locales son constantes para todos los elementos.
   %m_FFe = e_DatSet(iSet).m_FF;
   %siavare = e_DatSET(iSet).siavare;
   e_DatElemSet = e_DatSet(iSet).e_DatElem;
   eltype = e_DatElemSet.eltype;
   dofpe = e_DatElemSet.dofpe;
   npg = e_DatElemSet.npg;
   %npe = e_DatElemSet.npe;
   e_DatMatSet = e_DatSet(iSet).e_DatMat;
   %m_IndElemSet = e_DatSet(iSet).m_IndElemSet;
   m_NumElem = e_DatSet(iSet).m_NumElem;
   m_DefMacroSet = DefMacro{iSet}; %AA: La def macro es cte en todo el elemento aplicada en cada PG
   
   switch protype %AA
       case 0
            m_BT = e_DatSet(iSet).m_BT;
            m_DetJT = e_DatSet(iSet).m_DetJT;
%AA: movi.Estaba mas adelante. P/evitar mas if-elseif
%Como por convenci�n no deber�a usarse las news dentro de las funciones, solo cuando converja, no se las
%env�a como argumento, y por lo tanto se evita esa transferencia (tiene m�s influencia en el caso de
%cluster) de esas variables new. Por ello se debe generar en cada funci�n la porci�n que corresponda de
%las variables new (que tiene su costo, pero en el caso de un cluster es m�s conveniente), y tener en
%cuenta que si es multiescala es una estructura y si es material est�ndar es una matriz.
            sigma_old = e_VarEst_old(iSet).sigma;
            sigma_new   = e_VarEst_new(iSet).sigma;
            eps_new     = e_VarEst_new(iSet).eps;
            hvar_new    = e_VarEst_new(iSet).hvar;
            eps_fluct   = e_VarEst_new(iSet).eps_fluct;
            m_VarHistElemNew = e_VarEst_new(iSet).VarHistElem;
       case 1 %AA
            m_BT_d = e_DatSet(iSet).m_BT_d;
            m_DetJT_d = e_DatSet(iSet).m_DetJT_d;
            m_DerCa_p = e_DatSet(iSet).m_DerCa_p;
            m_DetJT_p = e_DatSet(iSet).m_DetJT_p;
            sigmaE_new   = e_VarEst_new(iSet).sigmaE; %Tensi�n efectiva
            sigmaT_new   = e_VarEst_new(iSet).sigmaT; %Tensi�n total
            sigmaE_old = e_VarEst_old(iSet).sigmaE; %AA Tensi�n efectiva
            sigmaT_old = e_VarEst_old(iSet).sigmaT; %AA Tensi�n efectiva
            eps_new     = e_VarEst_new(iSet).eps;
            hvar_new    = e_VarEst_new(iSet).hvar;
            eps_fluct   = e_VarEst_new(iSet).eps_fluct;
            m_VarHistElemNew = e_VarEst_new(iSet).VarHistElem;
   end %protype  %AA
   
   % Inicializaciones
   %aux_var = zeros(siavare,nElem);
   
   %La verificaci�n si esImplex es field de la estructura no ser�a necesario si a todos los modelos
   %constitutivos se agrega este campo.
   if isfield(e_DatMatSet,'esImplex')&&e_DatMatSet.esImplex
      %No se puede disminuir la matriz m_TensorTang cuando se calcula
      m_TensorTang = zeros(ntens,ntens,2*npg,nElem);
      %esImplex = 1;
   else
      m_TensorTang = zeros(ntens,ntens,npg,nElem);
      %esImplex = 0;
   end
   m_Ke = zeros(dofpe,dofpe,nElem);
   m_Fint = zeros(dofpe,nElem);
  
   % Tensor constitutivo el�stico
   %En este momento se est� precalculando en el read_data.
   %ce = c_elas(e_DatMatSet.young,e_DatMatSet.poiss,e_VG);
   
   % Coordenadas de los nodos de todos los elementos del set
   %Como se precalcula la matriz de deformaci�n, las coordenadas de los nodos no se necesita.
   %coord_n = reshape(f_CoordElem(xx,conec')',[],npe,nElem);
   coord_n = [];
   
   % Grados de libertad y coordenadas de los nodos de los elementos del set
   %dofElemSet = f_DofElem(reshape(conec',1,[]),ndn);
   dofElemSet = m_DofElem(:);
   %m_IndDof = repmat(dofElemSet,dofpe,1);
   %Ver si hay alguna forma de hacerlo m�s r�pido o mejor.
   %permute parece que hace copia de memoria, as� que es probable que sea mejor usar dos repmat.
   %m_Fil = reshape(permute(reshape(m_IndDof,dofpe,dofpe,[]),[2,1,3]),[],1);
   %m_Col = m_IndDof(:);
   m_FilFza = dofElemSet';
   m_Fil = reshape(repmat(reshape(dofElemSet,dofpe,[]),dofpe,1),1,[]);
   m_Col = reshape(repmat(m_FilFza,dofpe,1),1,[]);
   uElemSet  = reshape(u(dofElemSet),[],nElem);
   duElemSet = reshape(du(dofElemSet),[],nElem);
   %m_DefMacroSet = DefMacro(:,m_IndElemSet);
   
   %Solo para debug, se guarda el n�mero de set en e_VG (ya que esta variable se pasa a todas las
   %funciones).
   e_VG.iSet = iSet;
   
   %Se podr�a poner el parfor de dentro de cada eltype, as� evitar hacer la verificaci�n de tipo
   %de elemento para cada set.
   %parfor iElem = 1:nElem
   %for iElem = 1:nElem
   
   %Esta l�nea no puede estar si est� el parfor activado.
   %e_VG.iElem = iElem;
   %No se puede modificar una variable definida previa al loop, si no es sliced y si no es
   %interpretada como de reducci�n, ya que no sabe como interpretarla el MatLab (puede haber
   %superposici�n de resultados, al hacer reducci�n). Por eso cada Lab debe tener su "copia",
   %para modificarla en forma independiente.
   %Lo que puede ser lento es copiar para cada lab esa e_VG_Aux, que puede ser grande.
   %En realidad como cada procesador tiene su copia local del e_VG (ya que MatLab realiza una
   %copia por cada Lab de todas la variables, excepto las sliced, donde solo copia la parte
   %que le corresponde al procesador), y al "copiarse" esta al e_VG_Aux, lo que �nico que se
   %hace es copiarse el puntero, ya que no se est� modificando el e_VG_Aux. Luego lo que se
   %realiza es la modificaci�n de un valor de un campo (field), donde supuestamente
   %MatLab al cambiar un field de una estructura no hace copia de toda la estructura de nuevo,
   %si solo del field, por lo tanto no tendr�a que ser mucho m�s lenta.
   %Otra ser�a pasar la variable iElem como argumento en las funciones que llama dentro del
   %parfor.
   %e_VG_Aux = e_VG;
   %e_VG_Aux.iElem = m_IndElemSet(iElem);
   %e_VG_Aux.iElemSet = iElem;
   %fprintf('N�mero de elemento: %d\n',iElem)
   
   % Fuerza interna y tensor tangente del elemento
   switch protype %AA
       case 0
           switch eltype
              case 2
                 parfor iElem = 1:nElem
                 %for iElem = 1:nElem
                    e_VG_Aux2 = e_VG;
                    e_VG_Aux2.iElemSet = iElem;
                    e_VG_Aux2.iElemNum = m_NumElem(iElem);
                    %Tri�ngulo Est�ndar de 3 nodos
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem)] = ...
                       f_MatElem_tria_t1(...
                       coord_n,uElemSet(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,iElem),sigma_old(:,iElem),e_VG_Aux2);       % m_DefMacroSet(:,:,iElem) 
                 end
              case 4 %Cuadr�ngulo est�ndar de 4 nodos    
                 m_VolElem = e_DatSet(iSet).m_VolElem;
%                  parfor iElem = 1:nElem
                 for iElem = 1:nElem
                    e_VG_Aux4 = e_VG;
                    e_VG_Aux4.iElemSet = iElem;
                    e_VG_Aux4.iElemNum = m_NumElem(iElem);
                    %
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem)] =...
                       f_MatElem_quad_q1(...
                       uElemSet(:,iElem),eps_old(:,iElem), hvar_old(:,iElem),aux_var(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,iElem),sigma_old(:,iElem),...    %m_DefMacroSet(:,:,iElem)
                       m_VolElem(iElem),e_VG_Aux4);
                 end
              case 16    %AA: Cuadr�ngulo est�ndar de 8 nodos
                 m_VolElem = e_DatSet(iSet).m_VolElem;
%                  parfor iElem = 1:nElem
                 for iElem = 1:nElem
                    e_VG_Aux4 = e_VG;
                    e_VG_Aux4.iElemSet = iElem;
                    e_VG_Aux4.iElemNum = m_NumElem(iElem);
                    %
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem)] =...
                       f_MatElem_quad_q1(...
                       uElemSet(:,iElem),eps_old(:,iElem), hvar_old(:,iElem),aux_var(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,iElem),sigma_old(:,iElem),...    %m_DefMacroSet(:,:,iElem)
                       m_VolElem(iElem),e_VG_Aux4);
                 end
              case 31    % ELEMENTO TIPO BANDA : Cuadr�ngulo de 4 nodos
                 %               parfor iElem = 1:nElem
                 p_IntDissip  =  e_DatElemSet.pointersVHE.i_IntDissip  ;
                 p_IntEnergy  =  e_DatElemSet.pointersVHE.i_pIntEnergy;
                 m_IntDissip  =  m_VarHistElemOld (p_IntDissip ,:)  ;
                 m_pIntEnergy =  m_VarHistElemOld (p_IntEnergy ,:)  ;
                 m_pIntEnergy =  reshape(m_pIntEnergy, npg,[])  ;
                 ksb          =  e_DatElemSet.ksb;
                 %
                 parfor iElem = 1:nElem
                 %for iElem = 1:nElem
                    e_VG_Aux4 = e_VG;
                    e_VG_Aux4.iElemSet = iElem;
                    e_VG_Aux4.iElemNum = m_NumElem(iElem);
                    %
                    %m_incrDissip =  m_VarAuxElem (p_incrDissip   ,iElem)  ;
                    %
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem),m_IntDissip(iElem),...
                       m_pIntEnergy(:,iElem)] =...
                       f_MatElem_Banda_31(...
                       uElemSet(:,iElem),eps_old(:,iElem), hvar_old(:,iElem),aux_var(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,iElem),sigma_old(:,iElem),...    %m_DefMacroSet(:,:,iElem)
                       m_pIntEnergy(:,iElem), ksb(iElem),e_VG_Aux4);
                    %
                 end
                 m_VarHistElemNew (p_IntDissip ,: )      =  m_IntDissip(:)   ;
                 % m_pIntEnergy =  reshape(m_pIntEnergy,[],1)  ;
                 m_VarHistElemNew (p_IntEnergy ,: ) =  m_pIntEnergy   ;
              case 32
                 p_IntDissip  =  e_DatElemSet.pointersVHE.i_IntDissip  ;
                 p_IntEnergy  =  e_DatElemSet.pointersVHE.i_pIntEnergy;
                 m_IntDissip  =  m_VarHistElemOld (p_IntDissip ,:)  ;
                 m_pIntEnergy =  m_VarHistElemOld (p_IntEnergy ,:)  ;
                 m_pIntEnergy =  reshape(m_pIntEnergy, npg,[])  ;
                 ksb          =  e_DatElemSet.ksb;
                 parfor iElem = 1:nElem
        %        for iElem = 1:nElem
                    e_VG_Aux4 = e_VG;
                    e_VG_Aux4.iElemSet = iElem;
                    e_VG_Aux4.iElemNum = m_NumElem(iElem);
                    %
                    %m_incrDissip =  m_VarAuxElem (p_incrDissip   ,iElem)  ;
                    %
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem),m_IntDissip(iElem),...
                       m_pIntEnergy(:,iElem)] =...
                       f_MatElem_Banda_Tria1(...
                       uElemSet(:,iElem),eps_old(:,iElem), hvar_old(:,iElem),aux_var(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,iElem),sigma_old(:,iElem),...    %m_DefMacroSet(:,:,iElem)
                       m_pIntEnergy(:,iElem), ksb(iElem),e_VG_Aux4);
                    %
                  end
                  m_VarHistElemNew (p_IntDissip ,: )      =  m_IntDissip(:)   ;
                  %      m_pIntEnergy =  reshape(m_pIntEnergy,[],1)  ;
                  m_VarHistElemNew (p_IntEnergy ,: ) =  m_pIntEnergy   ;
              case 5
                 parfor iElem = 1:nElem
                    %                for iElem = 1:nElem
                    %Faltar�a juntar las funciones de las fuerzas internas y matriz de rigidez
                    %elemental. Tambi�n las funciones constitutivas.
                    %Falta precalcular las matrices de deformaciones y determinante de J. Esto est� hecho,
                    %falta pasarselo a las funciones.
                    [m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),hvar_new(:,iElem),...
                       aux_var(:,iElem)] = elem_int_force_barra2D(coord_n,uElemSet(:,iElem),...
                       hvar_old(:,iElem),e_DatMatEl,DefMacro_new,e_VG);
                    [m_Ke(:,:,iElem),m_TensorTang(:,:,:,iElem)] = elem_stiff_barra2D(coord_n,...
                       uElemSet(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),e_DatMatEl,ce,...
                       DefMacro_new,e_VG);
                 end
              case 7
                 error('Matrices Globales: Elemento finito no completo.')
              case 8      %Cuadr�ngulo BBar de 4 nodos
                 parfor iElem = 1:nElem
%                  for iElem = 1:nElem
                    e_VG_Aux8 = e_VG;
                    e_VG_Aux8.iElemSet = iElem;
                    e_VG_Aux8.iElemNum = m_NumElem(iElem);
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem)] = ...
                       f_MatElem_bbar_q1(...
                       uElemSet(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,iElem),sigma_old(:,iElem),e_VG_Aux8);%AA % m_DefMacroSet(:,:,iElem)
                 end
              case 10
                 m_resbT = zeros(ndime,nElem);
                 m_kbuT = zeros(ndime,dofpe,nElem);
                 m_invkbbT = zeros(ndime,ndime,nElem);
                 %Recuperaci�n de variable de salto (variable interna condesada)
                 m_Beta = c_GdlCond{iSet,1};
                 parfor iElem = 1:nElem
%                  for iElem = 1:nElem
                    e_VG_Aux10 = e_VG;
                    e_VG_Aux10.iElemSet = iElem;
                    e_VG_Aux10.iElemNum = m_NumElem(iElem);
                    %Tri�ngulo Est�ndar de 3 nodos con discontinuidad fuerte (SDA)
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),m_VarHistElemNew(:,iElem),...
                       aux_var(:,iElem),m_VarAuxElem(:,iElem),m_TensorTang(:,:,:,iElem),...
                       m_resbT(:,iElem),m_kbuT(:,:,iElem),m_invkbbT(:,:,iElem)] = ...
                       f_MatElem_SDA_tria_t1(...
                       uElemSet(:,iElem),m_Beta(:,iElem),hvar_old(:,iElem),m_VarHistElemOld(:,iElem),...
                       aux_var(:,iElem),m_VarAuxElem(:,iElem),e_DatElemSet,e_DatMatSet,...
                       m_BT(:,:,:,iElem),m_DetJT(:,iElem),m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),...
                       eps_old(:,iElem),e_VG_Aux10);
                 end
                 c_GdlCond{iSet,2} = m_resbT;
                 c_GdlCond{iSet,3} = m_kbuT;
                 c_GdlCond{iSet,4} = m_invkbbT;
              case 20    %Cuadr�ngulo de 4 nodos mixto con inyecci�n de deformaci�n
                 parfor iElem = 1:nElem
                    %                for iElem = 1:nElem
                    e_VG_Aux20 = e_VG;
                    e_VG_Aux20.iElemSet = iElem;
                    e_VG_Aux20.iElemNum = m_NumElem(iElem);
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),eps_fluct(:,iElem),...
                       hvar_new(:,iElem),m_VarHistElemNew(:,iElem),aux_var(:,iElem),m_TensorTang(:,:,:,iElem)] =...
                       f_MatElem_MixStrInj_quad_q1(...
                       uElemSet(:,iElem),hvar_old(:,iElem),m_VarHistElemOld(:,iElem),aux_var(:,iElem),...
                       m_VarAuxElem(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),e_VG_Aux20);
                 end
                 %##################################################################
              case 21    %Cuadrángulo de 4 nodos mixto con inyección de deformación
                 %##################################################################
                 m_Kbu_SDA      = zeros(ndime,dofpe,nElem);
                 m_KbbIN_SDA    = zeros(ndime,ndime,nElem);
                 m_Res_beta_SDA = zeros(ndime,nElem);
                 % MSL ELEMENT
                 m_Kbu_MSL      = zeros(ntens,dofpe,nElem);
                 m_KbbIN_MSL    = zeros(ntens,ntens,nElem);
                 m_Res_beta_MSL = zeros(ntens,nElem);
                 m_stressTilde_new = zeros(ntens*npg,nElem);

                 %Recuperaci�n de variable de salto (variable interna condesada)
                 Dbeta_SDA    = c_GdlCond{iSet,2};
                 Dbeta_MSL    = c_GdlCond{iSet,4};

                 i_indST        =  e_DatElemSet.pointersVHE.i_indST    ;
                 p_indActSTmacro=  e_DatElemSet.pointersVHE.p_indActSTmacro    ;
                 i_vectVHElem   =  e_DatElemSet.pointersVHE.i_vectVHElem ;
                 i_stressTilde  =  e_DatElemSet.pointersVHE.i_stressTilde ;

                 p_condBif      =  e_DatElemSet.pointersVAE.p_condBif ;
                 p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
                 p_leq_elem     =  e_DatElemSet.pointersVAE.p_leq_elem ;
                 p_kSD          =  e_DatElemSet.pointersVAE.p_kSD      ;
                 p_phi_grad     =  e_DatElemSet.pointersVAE.p_phi_grad  ;
                 p_n_tens       =  e_DatElemSet.pointersVAE.p_n_tens    ;
                 p_fii          =  e_DatElemSet.pointersVAE.p_fii       ;
                 p_injFactor    =  e_DatElemSet.pointersVAE.p_injFactor ;

                 NumsidesCutCPF = e_DatElemSet.NumsidesCutCPF;

                 m_indSTmacro_old  =  m_VarHistElemOld (i_indST        ,:)  ;
                 m_indActSTmacro_old =  m_VarHistElemOld (p_indActSTmacro        ,:)  ;
                 m_vectVHElem_old  =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                 m_indSTmacro      =  m_VarHistElemOld (i_indST        ,:)   ;
                 m_indActSTmacro     =  m_VarHistElemOld (p_indActSTmacro        ,:)  ;
                 m_vectVHElem      =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                 m_stressTilde_old =  m_VarHistElemOld (i_stressTilde   ,:)  ;         

                 elem_BifType  =    m_VarAuxElem(p_elem_BifType,:)  ;
                 kSD           =    m_VarAuxElem(p_kSD,:);

                 fact_inyect = e_DatElemSet.facIny;

                 parfor iElem = 1:nElem
%                  for iElem = 1:nElem

                    e_VG_Aux21 = e_VG;
                    e_VG_Aux21.iElemSet = iElem;
                    e_VG_Aux21.iElemNum = m_NumElem(iElem);


                    condBif         =  m_VarAuxElem(p_condBif,iElem) ;
                    leq_elem        =  m_VarAuxElem(p_leq_elem,iElem) ;
                    m_phi_grad      =  m_VarAuxElem(p_phi_grad,iElem) ;
                    m_phi_grad      =  reshape(m_phi_grad,4,2);
                    m_n_tens        =  m_VarAuxElem(p_n_tens,iElem) ;
                    m_n_tens        =  reshape(m_n_tens,4,2);

                    [m_Ke(:,:,iElem),m_Kbu_SDA(:,:,iElem),m_Kbu_MSL(:,:,iElem),m_KbbIN_SDA(:,:,iElem),...
                       m_KbbIN_MSL(:,:,iElem),m_Fint(:,iElem),m_Res_beta_SDA(:,iElem),...
                       m_Res_beta_MSL(:,iElem),sigma_new(:,iElem), ...
                       hvar_new(:,iElem),eps_new(:,iElem),m_TensorTang(:,:,:,iElem),...
                       m_indSTmacro(:,iElem) ,m_indActSTmacro(:,iElem),elem_BifType(iElem),...
                       kSD(iElem),m_vectVHElem(:,iElem),m_stressTilde_new(:,iElem)] = ...
                       f_MatElem_quad_q1_SDA ...
                       (duElemSet(:,iElem),eps_old(:,iElem),...
                       Dbeta_SDA(:,iElem),Dbeta_MSL(:,iElem),...
                       aux_var(:,iElem),condBif, leq_elem ,m_phi_grad,m_n_tens,sigma_old(:,iElem),...
                       hvar_old(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),...
                       m_DetJT(:,iElem),m_indSTmacro_old(iElem),m_indActSTmacro_old(iElem),...
                       e_VG_Aux21,m_vectVHElem_old(:,iElem),fact_inyect,...
                       NumsidesCutCPF(iElem),m_stressTilde_old(:,iElem));

                 end

                 m_VarHistElemNew (i_indST   ,:)      = m_indSTmacro  ;
                 m_VarHistElemNew (p_indActSTmacro   ,:) = m_indActSTmacro  ;
                 m_VarHistElemNew (i_vectVHElem   ,:) = m_vectVHElem  ;
                 m_VarHistElemNew (i_stressTilde  ,:) = m_stressTilde_new  ;

                 m_VarAuxElem(p_elem_BifType,:)       =  elem_BifType  ;
                 m_VarAuxElem(p_kSD,:)                =  kSD  ;



                 % SDA ELEMENT
                 c_GdlCond{iSet,5}  =   m_Kbu_SDA ;
                 c_GdlCond{iSet,6}  =   m_KbbIN_SDA ;
                 c_GdlCond{iSet,7}  =   m_Res_beta_SDA;
                 % MSL ELEMENT
                 c_GdlCond{iSet,8}  =   m_Kbu_MSL;
                 c_GdlCond{iSet,9}  =   m_KbbIN_MSL ;
                 c_GdlCond{iSet,10} =   m_Res_beta_MSL;

                 %##################################################################
              case 22    %Cuadrangulo de 4 nodos mixto con inyeccion de deformacion
                 % MODELO SANTA FE
                 %##################################################################
                 m_Kbu_SDA      = zeros(ndime,dofpe,nElem);
                 m_KbbIN_SDA    = zeros(ndime,ndime,nElem);
                 m_Res_beta_SDA = zeros(ndime,nElem);
                 % MSL ELEMENT
                 m_Kbu_MSL      = zeros(ntens,dofpe,nElem);
                 m_KbbIN_MSL    = zeros(ntens,ntens,nElem);
                 m_Res_beta_MSL = zeros(ntens,nElem);
                 m_stressTilde_new = zeros(ntens*npg,nElem);

                 %Recuperaci�n de variable de salto (variable interna condesada)
                 Dbeta_SDA    = c_GdlCond{iSet,2};
                 Dbeta_MSL    = c_GdlCond{iSet,4};

                 i_indST        =  e_DatElemSet.pointersVHE.i_indST    ;
                 p_indActSTmacro=  e_DatElemSet.pointersVHE.p_indActSTmacro    ;
                 i_vectVHElem   =  e_DatElemSet.pointersVHE.i_vectVHElem ;
                 i_stressTilde  =  e_DatElemSet.pointersVHE.i_stressTilde ;

                 p_condBif      =  e_DatElemSet.pointersVAE.p_condBif ;
                 p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
                 p_leq_elem     =  e_DatElemSet.pointersVAE.p_leq_elem ;
                 p_kSD          =  e_DatElemSet.pointersVAE.p_kSD      ;
                 p_phi_grad     =  e_DatElemSet.pointersVAE.p_phi_grad  ;
                 p_n_tens       =  e_DatElemSet.pointersVAE.p_n_tens    ;

                 NumsidesCutCPF = e_DatElemSet.NumsidesCutCPF;

                 m_indSTmacro_old  =  m_VarHistElemOld (i_indST        ,:)  ;
                 m_indActSTmacro_old =  m_VarHistElemOld (p_indActSTmacro        ,:)  ;
                 m_vectVHElem_old  =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                 m_indSTmacro      =  m_VarHistElemOld (i_indST        ,:)   ;
                 m_indActSTmacro     =  m_VarHistElemOld (p_indActSTmacro        ,:)  ;
                 m_vectVHElem      =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                 m_stressTilde_old =  m_VarHistElemOld (i_stressTilde   ,:)  ;

                 elem_BifType  =    m_VarAuxElem(p_elem_BifType,:)  ;
                 kSD           =    m_VarAuxElem(p_kSD,:);

                 fact_inyect = e_DatElemSet.facIny;

                 parfor iElem = 1:nElem
%                  for iElem = 1:nElem

                    e_VG_Aux21 = e_VG;
                    e_VG_Aux21.iElemSet = iElem;
                    e_VG_Aux21.iElemNum = m_NumElem(iElem);


                    condBif         =  m_VarAuxElem(p_condBif,iElem) ;
                    leq_elem        =  m_VarAuxElem(p_leq_elem,iElem) ;
                    m_phi_grad      =  m_VarAuxElem(p_phi_grad,iElem) ;
                    m_phi_grad      =  reshape(m_phi_grad,4,2);
                    m_n_tens        =  m_VarAuxElem(p_n_tens,iElem) ;
                    m_n_tens        =  reshape(m_n_tens,4,2);

                    [m_Ke(:,:,iElem),m_Kbu_SDA(:,:,iElem),m_Kbu_MSL(:,:,iElem),m_KbbIN_SDA(:,:,iElem),...
                       m_KbbIN_MSL(:,:,iElem),m_Fint(:,iElem),m_Res_beta_SDA(:,iElem),...
                       m_Res_beta_MSL(:,iElem),sigma_new(:,iElem), ...
                       hvar_new(:,iElem),eps_new(:,iElem),m_TensorTang(:,:,:,iElem),...
                       m_indSTmacro(:,iElem) ,m_indActSTmacro(:,iElem),elem_BifType(iElem),...
                       kSD(iElem),m_vectVHElem(:,iElem),m_stressTilde_new(:,iElem)] = ...
                       f_MatElem_quad_q1_SDA_STAFE(...
                       duElemSet(:,iElem),eps_old(:,iElem),...
                       Dbeta_SDA(:,iElem),Dbeta_MSL(:,iElem),...
                       aux_var(:,iElem),condBif, leq_elem ,m_phi_grad,m_n_tens,sigma_old(:,iElem),...
                       hvar_old(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),...
                       m_DetJT(:,iElem),m_indSTmacro_old(iElem),m_indActSTmacro_old(iElem),...
                       e_VG_Aux21,m_vectVHElem_old(:,iElem),fact_inyect,...
                       NumsidesCutCPF(iElem),m_stressTilde_old(:,iElem));

                 end

                 m_VarHistElemNew (i_indST   ,:)      = m_indSTmacro  ;
                 m_VarHistElemNew (p_indActSTmacro   ,:) = m_indActSTmacro  ;
                 m_VarHistElemNew (i_vectVHElem   ,:) = m_vectVHElem  ;
                 m_VarHistElemNew (i_stressTilde  ,:) = m_stressTilde_new  ;

                 m_VarAuxElem(p_elem_BifType,:)       =  elem_BifType  ;
                 m_VarAuxElem(p_kSD,:)                =  kSD  ;


                 % SDA ELEMENT
                 c_GdlCond{iSet,5}  =   m_Kbu_SDA ;
                 c_GdlCond{iSet,6}  =   m_KbbIN_SDA ;
                 c_GdlCond{iSet,7}  =   m_Res_beta_SDA;
                 % MSL ELEMENT
                 c_GdlCond{iSet,8}  =   m_Kbu_MSL;
                 c_GdlCond{iSet,9}  =   m_KbbIN_MSL ;
                 c_GdlCond{iSet,10} =   m_Res_beta_MSL;
                 %##################################################################
              case 23    %Cuadrangulo de 4 nodos mixto SD
                 % MODELO SANTA FE
                 %##################################################################

                 m_VolElem =  e_DatSet(iSet).m_VolElem;

                 m_Kbu_SDA      = zeros(ndime,dofpe,nElem);
                 m_KbbIN_SDA    = zeros(ndime,ndime,nElem);
                 m_Res_beta_SDA = zeros(ndime,nElem);
                 % MSL ELEMENT
                 m_Kbu_MSL      = zeros(ntens,dofpe,nElem);
                 m_KbbIN_MSL    = zeros(ntens,ntens,nElem);
                 m_Res_beta_MSL = zeros(ntens,nElem);
                 m_stressTilde_new = zeros(ntens*npg,nElem);

                 %Recuperación de variable de salto (variable interna condesada)
                 Dbeta_SDA    = c_GdlCond{iSet,2};
                 Dbeta_MSL    = c_GdlCond{iSet,4};

                 i_indST        =  e_DatElemSet.pointersVHE.i_indST    ;
                 p_indActSTmacro=  e_DatElemSet.pointersVHE.p_indActSTmacro    ;
                 i_vectVHElem   =  e_DatElemSet.pointersVHE.i_vectVHElem ;
                 i_stressTilde  =  e_DatElemSet.pointersVHE.i_stressTilde ;

                 p_condBif      =  e_DatElemSet.pointersVAE.p_condBif ;
                 p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
                 p_leq_elem     =  e_DatElemSet.pointersVAE.p_leq_elem ;
                 p_kSD          =  e_DatElemSet.pointersVAE.p_kSD      ;
                 p_phi_grad     =  e_DatElemSet.pointersVAE.p_phi_grad  ;
                 p_n_tens       =  e_DatElemSet.pointersVAE.p_n_tens    ;
                 p_CeBif        =  e_DatElemSet.pointersVAE.pCeBif;
                 p_Traccion = e_DatElemSet.pointersVAE.p_Traccion;

                 NumsidesCutCPF = e_DatElemSet.NumsidesCutCPF;

                 m_indSTmacro_old  =  m_VarHistElemOld (i_indST        ,:)  ;
                 m_indActSTmacro_old =  m_VarHistElemOld(p_indActSTmacro,:)  ;
                 m_vectVHElem_old  =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                 m_indSTmacro      =  m_VarHistElemOld (i_indST        ,:)  ;
                 m_indActSTmacro   =  m_VarHistElemOld (p_indActSTmacro,:)  ;
                 m_vectVHElem      =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                 m_stressTilde_old =  m_VarHistElemOld (i_stressTilde   ,:)  ;

                 elem_BifType  =    m_VarAuxElem(p_elem_BifType,:)  ;
                 kSD           =    m_VarAuxElem(p_kSD,:);
                 m_CeBif = reshape(m_VarAuxElem(p_CeBif,:),ntens,ntens,[]);
                 %En realidad no se env�a como argumento de entrada, pero como se necesita almacenar la salida para
                 %cada elemento, se la inicializa.
                 m_Tracc = zeros(ndime,nElem);

                 %o_Par = Par(e_VG.nLab);
                 parfor iElem = 1:nElem
%                  for iElem = 1:nElem
                    %Par.tic;

                    e_VG_Aux23 = e_VG;
                    e_VG_Aux23.iElemSet = iElem;
                    e_VG_Aux23.iElemNum = m_NumElem(iElem);            

                    %Se crea esta variable adicional para que el MatLab se d� cuenta que m_VarAuxElem es una
                    %variable slice. Ver que, por ejemplo, en m_VarAuxElem(p_condBif,iElem), como no sabe como indexa
                    %el vector p_condBif dentro de m_VarAuxElem no puede definir que es slice.
                    m_VarAuxiElem = m_VarAuxElem(:,iElem);
                    condBif         =  m_VarAuxiElem(p_condBif) ;
                    leq_elem        =  m_VarAuxiElem(p_leq_elem) ;
                    m_phi_grad      =  m_VarAuxiElem(p_phi_grad) ;
                    m_phi_grad      =  reshape(m_phi_grad,ntens,2);
                    m_n_tens        =  m_VarAuxiElem(p_n_tens) ;
                    m_n_tens        =  reshape(m_n_tens,ntens,2);

                    [m_Ke(:,:,iElem),m_Kbu_SDA(:,:,iElem),m_Kbu_MSL(:,:,iElem),m_KbbIN_SDA(:,:,iElem),...
                       m_KbbIN_MSL(:,:,iElem),m_Fint(:,iElem),m_Res_beta_SDA(:,iElem),...
                       m_Res_beta_MSL(:,iElem),sigma_new(:,iElem), ...
                       hvar_new(:,iElem),eps_new(:,iElem),m_TensorTang(:,:,:,iElem),m_CeBif(:,:,iElem),...
                       m_indSTmacro(:,iElem) ,m_indActSTmacro(:,iElem),...
                       kSD(iElem),m_vectVHElem(:,iElem),m_stressTilde_new(:,iElem),m_Tracc(:,iElem)] = ...
                       f_MatElem_quad_q23_SDA_STAFE(...
                       duElemSet(:,iElem),eps_old(:,iElem),...
                       Dbeta_SDA(:,iElem),Dbeta_MSL(:,iElem),...
                       aux_var(:,iElem),condBif, leq_elem ,m_phi_grad,m_n_tens,m_CeBif(:,:,iElem),...
                       sigma_old(:,iElem),...
                       hvar_old(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),...
                       m_DetJT(:,iElem),m_indSTmacro_old(iElem),m_indActSTmacro_old(iElem),...
                       e_VG_Aux23,m_vectVHElem_old(:,iElem),...
                       NumsidesCutCPF(iElem),m_stressTilde_old(:,iElem),m_VolElem(iElem));

                    %o_Par(iElem) = Par.toc;            
                 end
                 %stop(o_Par);

                 m_VarHistElemNew (i_indST   ,:)      = m_indSTmacro  ;
                 m_VarHistElemNew (p_indActSTmacro   ,:) = m_indActSTmacro  ;
                 m_VarHistElemNew (i_vectVHElem   ,:) = m_vectVHElem  ;
                 m_VarHistElemNew (i_stressTilde  ,:) = m_stressTilde_new  ;

                 m_VarAuxElem(p_elem_BifType,:)       =  elem_BifType  ;
                 m_VarAuxElem(p_kSD,:)                =  kSD  ;
                 m_VarAuxElem(p_CeBif,:) = reshape(m_CeBif,ntens*ntens,[]);
                 %Se asume que en m_Tracc viene las tracciones totales
                 m_VarAuxElem(p_Traccion,:) = m_Tracc;


                 % SDA ELEMENT
                 c_GdlCond{iSet,5}  =   m_Kbu_SDA ;
                 c_GdlCond{iSet,6}  =   m_KbbIN_SDA ;
                 c_GdlCond{iSet,7}  =   m_Res_beta_SDA;
                 % MSL ELEMENT
                 c_GdlCond{iSet,8}  =   m_Kbu_MSL;
                 c_GdlCond{iSet,9}  =   m_KbbIN_MSL ;
                 c_GdlCond{iSet,10} =   m_Res_beta_MSL;

              %
              %LARGE DEFORMATIONS
              case 108  %Cuadr�ngulo de 4 nodos FBar.
                 %m_VolElem = e_DatSet(iSet).m_VolElem;         
                 p_FreeEnergy  =  1+[1:npg*4] ; %e_DatElemSet.pointersVHE.i_pIntEnergy;   % +3 para llevar Psi_e_vol, Psi_e_dev, Psi_p
                 m_FreeEnergy =  m_VarHistElemOld (p_FreeEnergy ,:)  ;
                 m_FreeEnergy =  reshape(m_FreeEnergy, npg*4,[])  ;

                 parfor iElem = 1:nElem %JLM
%                   for iElem = 1:nElem
                    e_VG_Aux108 = e_VG;
                    e_VG_Aux108.iElemSet = iElem;
                    e_VG_Aux108.iElemNum = m_NumElem(iElem);
                    %
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem),m_FreeEnergy(:,iElem)] =...
                       f_MatElem_FBar_q1(...
                       uElemSet(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),m_VarAuxElem(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,iElem),eps_old(:,iElem),sigma_old(:,iElem),...
                       e_VG_Aux108);
                    %
                    % C�lculo del Ke num�rico
        %             duAux108 = 1e-8;
        %             uAux108 = uElemSet(:,iElem);
        %             m_udf = uAux108;
        %             kt = zeros(dofpe,dofpe);
        %             for i = 1:dofpe
        %                m_udf(i) = uAux108(i)+duAux108;
        %                [~,F1] = f_MatElem_FBar_q1(...
        %                   m_udf,hvar_old(:,iElem),aux_var(:,iElem),m_VarAuxElem(:,iElem),...
        %                   e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
        %                   m_DefMacroSet(:,iElem),eps_old(:,iElem),sigma_old(:,iElem),...
        %                   e_VG_Aux108);
        %                m_udf(i) = uAux108(i)-duAux108;
        %                [~,F0] = f_MatElem_FBar_q1(...
        %                   m_udf,hvar_old(:,iElem),aux_var(:,iElem),m_VarAuxElem(:,iElem),...
        %                   e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
        %                   m_DefMacroSet(:,iElem),eps_old(:,iElem),sigma_old(:,iElem),...
        %                   e_VG_Aux108);
        %                kt(:,i) = (F1-F0)/(2*duAux108);
        %                m_udf(i) = uAux108(i);
        %             end
        % %             %norm(m_Ke(:,:,iElem)-kt)/norm(kt)
        %             m_Ke(:,:,iElem) = kt;

        %             m_VarHistElemNew(p_FreeEnergy ,: ) =  m_FreeEnergy ;
                 end
                 m_VarHistElemNew(p_FreeEnergy ,: ) =  m_FreeEnergy ;

           end %eltype %AA
           
       case 1 %AA: add toda la construccion
           switch eltype
               case 16    %AA: Cuadr�ngulo est�ndar de 8 nodos (Aplicado p/medio bifase)
                 %m_VolElem_d = e_DatSet(iSet).m_VolElem_d;
                 %m_VolElem_p = e_DatSet(iSet).m_VolElem_p;
                 m_FF_p = e_DatSet(iSet).m_FF_p;
%                  parfor iElem = 1:nElem
                 for iElem = 1:nElem
                    e_VG_Aux16 = e_VG;
                    e_VG_Aux16.iElemSet = iElem;
                    e_VG_Aux16.iElemNum = m_NumElem(iElem);
                    %
                    [m_Ke(:,:,iElem),m_Fint(:,iElem),sigmaE_new(:,iElem),sigmaT_new(:,iElem),...
                       eps_new(:,iElem),eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem)] =...
                       f_MatElem_Bifase(...
                       uElemSet(:,iElem),eps_old(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT_d(:,:,:,iElem),m_DetJT_d(:,iElem),...
                       m_DerCa_p(:,:,:,iElem),m_DetJT_p(:,iElem),m_FF_p,m_DefMacroSet(:,iElem),...  %m_DefMacroSet(:,:,iElem)
                       sigmaE_old(:,iElem),sigmaT_old(:,iElem),e_VG_Aux16);
                 end %for(iElem)
               otherwise %eltype
                     error('Elemento no desarrollado para medio bifasico')
           end %eltype %AA
           
   end %protype %AA
   
   if protype==0 %AA
           e_VarEst_new(iSet).sigma = sigma_new;
           e_VarEst_new(iSet).eps = eps_new;
           e_VarEst_new(iSet).hvar = hvar_new;
           e_VarEst_new(iSet).eps_fluct = eps_fluct;
           e_VarEst_new(iSet).VarHistElem = m_VarHistElemNew;
           e_VarAux(iSet).VarAuxGP = aux_var;
           e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
           c_Ke{iSet} = m_Ke(:);
           c_Fint{iSet} = m_Fint(:);
           c_Fil{iSet} = m_Fil;
           c_Col{iSet} = m_Col;
           c_FilFza{iSet} = m_FilFza;
           c_TensorTang{iSet} = m_TensorTang;
   elseif protype==1
           e_VarEst_new(iSet).sigmaE = sigmaE_new;
           e_VarEst_new(iSet).sigmaT = sigmaT_new;
           e_VarEst_new(iSet).eps = eps_new;
           e_VarEst_new(iSet).hvar = hvar_new;
           e_VarEst_new(iSet).eps_fluct = eps_fluct;
           e_VarEst_new(iSet).VarHistElem = m_VarHistElemNew;
           e_VarAux(iSet).VarAuxGP = aux_var;
           e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
           c_Ke{iSet} = m_Ke(:);
           c_Fint{iSet} = m_Fint(:);
           c_Fil{iSet} = m_Fil;
           c_Col{iSet} = m_Col;
           c_FilFza{iSet} = m_FilFza;
           c_TensorTang{iSet} = m_TensorTang;
   end %if(protype) %AA
   
end %for(iSet)

% Ensamble de matriz de fuerzas internas global
Fint = sparse([c_FilFza{:}],1,cat(1,c_Fint{:}));

% Ensamble de matriz de rigidez global
KT = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:}));

end
