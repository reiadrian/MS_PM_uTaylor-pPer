%Esta funcion permite calcular la matriz de rigidez y el vector de fuerzas
%internas de cada elemento y global
function [KT,Fint,c_GdlCond,e_VarEst_new,e_VarAux,c_sig_eps,o_Par] = ...
    f_Par_MatGlobalesPM(xx,u,u_old,c_GdlCond,e_DatSet,e_VarEst_new,...
    e_VarEst_old,e_VarAux,DefMacro,GradPorMacro_dup,...
    GradPorMacro_up_n,PorMacro,DefMacro_new,GradPorMacro_up_new,PorMacro_new,e_VG)

o_Par = [];
%Para no transferir estas variables que no se utilizan dentro del parfor.
e_VG.smooth_alpha = [];
e_VG.smooth_dalpha = [];
e_VG.smooth_dalpha_old = [];
e_VG.MGlobInv = [];

% Recupera variables globales
ntens = e_VG.ntens;
nSet = e_VG.nSet;

%############################################################################################
% INICIALIZACIONES
%Matriz tangente
c_Ke = cell(nSet,1);
%Vector de fuerzas internas
c_Fint = cell(nSet,1);
%Vector que permite armar el vector la matriz tangente como sparse
c_Fil = cell(nSet,1);
%Vector que permite armar el vector la matriz tangente como sparse
c_Col = cell(nSet,1);
%Vector que permite armar el vector de fuerzas internas como sparse
c_FilFza = cell(nSet,1);
% c_TensorTang = cell(nSet,1);


%OPERADORES TANGENTES
%Operador tangente constitutivo dsigma/deps
c_sig_eps = cell(nSet,1);
%Operador tangente constitutivo dchi/deps
c_chi_eps= cell(nSet,1);
%Operador tangente constitutivo dV/deps
c_V_eps= cell(nSet,1);
%Operador tangenteconstitutivo dsigma/dp
c_sig_p= cell(nSet,1);
%Operador tangenteconstitutivo dsigma/dp
c_chi_p= cell(nSet,1);
%Operador tangente constitutivo dchi/dp  
c_V_p= cell(nSet,1);
%Operador tangente constitutivo dsigma/dphi 
c_sig_phi= cell(nSet,1);
%Operador tangente constitutivo dchi/dphi 
c_chi_phi= cell(nSet,1);
%Operador tangente constitutivo dV/dphi   
c_V_phi= cell(nSet,1);
%############################################################################################

%Loop sobre el set de elementos (tipo y material del elemento)
for iSet = 1:nSet
    % Recuperacion de variables
    %(evita llamar desde las estructura de los sets en el bucle de los elementos, ver que solo
    %recupera punteros)
    %Numero total de elementos
    nElem = e_DatSet(iSet).nElem;
    %Matriz con grados de libertad
    m_DofElem = e_DatSet(iSet).m_DofElem;
    
    %Variables historicas de cada modelo constitutivo al tiempo "n"
    hvar_old = e_VarEst_old(iSet).hvar;
    aux_var = e_VarAux(iSet).VarAuxGP;
    m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
    %Estructura de datos del elemento finito aplicado
    e_DatElemSet = e_DatSet(iSet).e_DatElem;
    %Estructura de datos del material
    e_DatMatSet = e_DatSet(iSet).e_DatMat;
    %Grados de libertad conforme al elemento adoptado
    dofpe = e_DatElemSet.dofpe;
    %Numero de puntos de Gauss
    npg = e_DatElemSet.npg;
    
    m_NumElem = e_DatSet(iSet).m_NumElem;
    
    %Delta de macro-deformaciones
    m_DefMacroSet = DefMacro{iSet}; 
    %Delta de macro-gradiente de poro presiones 
    m_GradPorMacroSet_dup = GradPorMacro_dup{iSet}; 
    %Macro-gradiente de poro presiones al tiempo "n"
    m_GradPorMacroSet_up_n = GradPorMacro_up_n{iSet};
    %Delta de macro-poro presiones
    m_PorMacroSet = PorMacro{iSet}; 
    
    %Macro-deformaciones al tiempo "n+1"
    m_DefMacroSet_new = DefMacro_new{iSet}; 
    %Macro-gradiente de poro presiones al tiempo "n+1"
    m_GradPorMacroSet_up_new = GradPorMacro_up_new{iSet};
    %Macro-poro presiones  al tiempo "n+1"
    m_PorMacroSet_new = PorMacro_new{iSet}; 
    
    %Matriz B de desplazamientos-deformaciones
    m_BT_d = e_DatSet(iSet).m_BT_d; 
    %Determinante del Jacobiano en desplazamientos
    m_DetJT_d = e_DatSet(iSet).m_DetJT_d; 
    %Matriz en derivadas cartesianas para poro presiones
    m_DerCa_p = e_DatSet(iSet).m_DerCa_p; 
    %Determinante del Jacobiano en poro presiones
    m_DetJT_p = e_DatSet(iSet).m_DetJT_p; 
    %Funciones de interpolacion del elemento cuadrilatero bicuadratico
    %(desplazamientos)
    m_FF_d = e_DatSet(iSet).m_FF_d;
    %Funciones de interpolacion del elemento cuadrilatero bilineal (poro
    %presiones)
    m_FF_p = e_DatSet(iSet).m_FF_p;
    %Matriz de conectividades
    conec = e_DatSet(iSet).conec;

    %EN DESUSO AUN
    eps_old = e_VarEst_old(iSet).eps;
    sigmaE_old = e_VarEst_old(iSet).sigmaE; 
    sigmaT_old = e_VarEst_old(iSet).sigmaT; 

    %############################################################################################
    %VARIABLES PRIMITIVAS
    %Delta de deformaciones
    eps_new= e_VarEst_new(iSet).eps;
    %Delta de gradiente de poro presiones
    phi_new= e_VarEst_new(iSet).phi;
    %Delta de poro presiones
    porpr_new= e_VarEst_new(iSet).porpr;
    
    %VARIABLES PRIMITIVAS FLUCTUANTES
    %Delta de micro-deformaciones fluctuantes o delta de macro-deformaciones
    eps_fluct= e_VarEst_new(iSet).eps_fluct;
    %Delta de gradiente de micro-poro presiones fluctuantes o delta de gradiente de macro-poro presiones
    phi_fluct= e_VarEst_new(iSet).phi_fluct;
    %Delta de micro-poro presiones fluctuantes o delta de macro-poro presiones
    p_fluct= e_VarEst_new(iSet).p_fluct;
    
    %VARIABLES DUALES
    %Delta de tensiones efectivas
    sigmaE_new= e_VarEst_new(iSet).sigmaE;
    %Delta de tensiones totales
    sigmaT_new= e_VarEst_new(iSet).sigmaT;
    %Delta de tasa del contenido de masa del fluido
    mflu_new= e_VarEst_new(iSet).mflu;

    %Componente estatica de la velocidad de filtracion  al tiempo "n+theta"
    velflu_sta= e_VarEst_new(iSet).velflu_sta;
    %Delta de la componente dinamica de la velocidad de filtracion
    velflu_total= e_VarEst_new(iSet).velflu_total;
    
    %Velocidad de filtracion  al tiempo "n+theta"
    velflu_new= e_VarEst_new(iSet).velflu;
    
    %Variables historicas de cada modelo constitutivo al tiempo "n+1"    
    hvar_new    = e_VarEst_new(iSet).hvar;           
    m_VarHistElemNew = e_VarEst_new(iSet).VarHistElem;

    % Inicializaciones
    %La verificacion si esImplex es field de la estructura no seria necesario si a todos los modelos
    %constitutivos se agrega este campo.
    if isfield(e_DatMatSet,'esImplex') && e_DatMatSet.esImplex 
        %No se puede disminuir la matriz m_TensorTang cuando se calcula 
        %############################################################################################
        %OPERADORES TANGENTES
        %Operador tangente constitutivo dsigma/deps
        m_sig_eps = zeros(ntens,ntens,2*npg,nElem);
        %esImplex = 1;
    else
        %############################################################################################
        %OPERADORES TANGENTES
        %Operador tangente constitutivo dsigma/deps
        m_sig_eps= zeros(ntens,ntens,npg,nElem);
        %Operador tangente constitutivo dchi/deps
        m_chi_eps= zeros(1,ntens,npg,nElem);
        %Operador tangente constitutivo dV/deps
        m_V_eps= zeros(2,ntens,npg,nElem);
        %Operador tangenteconstitutivo dsigma/dp
        m_sig_p= zeros(ntens,1,npg,nElem);
        %Operador tangente constitutivo dchi/dp 
        m_chi_p= zeros(1,1,npg,nElem);
        %Operador tangente constitutivo dV/dp
        m_V_p=zeros(2,1,npg,nElem);
        %Operador tangente constitutivo dsigma/dphi
        m_sig_phi= zeros(ntens,2,npg,nElem);
        %Operador tangente constitutivo dchi/dphi   
        m_chi_phi= zeros(1,2,npg,nElem);
        %Operador tangente constitutivo dV/dphi 
        m_V_phi= zeros(2,2,npg,nElem);
        %esImplex = 0;
    end
     %############################################################################################
    if e_VG.conshyp==15 %Modelo micro-escala con restricciones de media volumetrica 
        %nula micro-poro presiones fluctuantes
        %Numero de grados de libertad totales
        ndoft=e_VG.ndoft;
        %Por el grado de libertad adicional del multiplicador de Lagrange
        %CONSTANTE. Agrego solo aca porque entiendo requiere menos
        %modificaciones que si lo hiciera en read_data
        %Matriz de rigidez del elemento con filas y columnas adicionales
        %para captar las restricciones
        m_Ke = zeros(dofpe+1,dofpe+1,nElem);
        %Vector de fuerzas internas del elemento con filas adicionales
        %para captar las restricciones
        m_Fint = zeros(dofpe+1,nElem); 
        
        % Grados de libertad y coordenadas de los nodos de los elementos del set
        dofElemSet_p = [m_DofElem; repmat(ndoft,1,nElem)];
        dofElemSet = dofElemSet_p(:);
        m_FilFza = dofElemSet';
        m_Fil = reshape(repmat(reshape(dofElemSet,dofpe+1,[]),dofpe+1,1),1,[]);
        m_Col = reshape(repmat(m_FilFza,dofpe+1,1),1,[]);
        %Vector de incognitas al tiempo  "n+1" (Desplazamientos, poro
        %presiones y multiplicadores de Lagrange)
        uElemSet  = reshape(u(dofElemSet),[],nElem);% repmat(u(ndoft),1,nElem)];
        %Vector de incognitas al tiempo  "n" (Desplazamientos, poro
        %presiones y multiplicadores de Lagrange)
        u_oldElemSet = reshape(u_old(dofElemSet),[],nElem);% repmat(u(ndoft),1,nElem)];
    
    elseif e_VG.conshyp==16  %Modelo micro-escala con restricciones de media volumetrica 
        %nula en micro-desplazamientos fluctuantes 
        %Numero de grados de libertad totales
        ndoft=e_VG.ndoft;
        %Por el grado de libertad adicional del multiplicador de Lagrange
        %CONSTANTE. Agrego solo aca porque entiendo requiere menos
        %modificaciones que si lo hiciera en read_data
        %Matriz de rigidez del elemento con filas y columnas adicionales
        %para captar las restricciones
        m_Ke = zeros(dofpe+2,dofpe+2,nElem);
        %Vector de fuerzas internas del elemento con filas adicionales
        %para captar las restricciones
        m_Fint = zeros(dofpe+2,nElem); 
        
        % Grados de libertad y coordenadas de los nodos de los elementos del set
        dofElemSet_p = [m_DofElem; repmat(ndoft-1,1,nElem); repmat(ndoft,1,nElem)];
        dofElemSet = dofElemSet_p(:);
        m_FilFza = dofElemSet';
        m_Fil = reshape(repmat(reshape(dofElemSet,dofpe+2,[]),dofpe+2,1),1,[]);
        m_Col = reshape(repmat(m_FilFza,dofpe+2,1),1,[]);
        %Vector de incognitas al tiempo  "n+1" (Desplazamientos, poro
        %presiones y multiplicadores de Lagrange)
        uElemSet  = reshape(u(dofElemSet),[],nElem);% repmat(u(ndoft),1,nElem)];
        %Vector de incognitas al tiempo  "n" (Desplazamientos, poro
        %presiones y multiplicadores de Lagrange)
        u_oldElemSet = reshape(u_old(dofElemSet),[],nElem);% repmat(u(ndoft),1,nElem)];
     elseif e_VG.conshyp==17 %Modelo micro-escala con restricciones de media volumetrica 
        %nula en micro-desplazamientos fluctuantes y micro-poro presiones fluctuantes
        %Numero de grados de libertad totales
        ndoft=e_VG.ndoft;
        %Por el grado de libertad adicional del multiplicador de Lagrange
        %CONSTANTE. Agrego solo aca porque entiendo requiere menos
        %modificaciones que si lo hiciera en read_data
        %Matriz de rigidez del elemento con filas y columnas adicionales
        %para captar las restricciones
        m_Ke = zeros(dofpe+3,dofpe+3,nElem);
        %Vector de fuerzas internas del elemento con filas adicionales
        %para captar las restricciones
        m_Fint = zeros(dofpe+3,nElem); 
        
        % Grados de libertad y coordenadas de los nodos de los elementos del set
        dofElemSet_p = [m_DofElem; repmat(ndoft-2,1,nElem); repmat(ndoft-1,1,nElem); repmat(ndoft,1,nElem)];
        dofElemSet = dofElemSet_p(:);
        m_FilFza = dofElemSet';
        m_Fil = reshape(repmat(reshape(dofElemSet,dofpe+3,[]),dofpe+3,1),1,[]);
        m_Col = reshape(repmat(m_FilFza,dofpe+3,1),1,[]);
        %Vector de incognitas al tiempo  "n+1" (Desplazamientos, poro
        %presiones y multiplicadores de Lagrange)
        uElemSet  = reshape(u(dofElemSet),[],nElem);% repmat(u(ndoft),1,nElem)];
        %Vector de incognitas al tiempo  "n" (Desplazamientos, poro
        %presiones y multiplicadores de Lagrange)
        u_oldElemSet = reshape(u_old(dofElemSet),[],nElem);% repmat(u(ndoft),1,nElem)];
    
    else %Cualquier otro modelo constitutivo macro o micro-escala
        %Matriz de rigidez del elemento
        m_Ke = zeros(dofpe,dofpe,nElem);
        %Vector de fuerzas internas del elemento
        m_Fint = zeros(dofpe,nElem); 
        
        % Grados de libertad y coordenadas de los nodos de los elementos del set
        dofElemSet = m_DofElem(:);
        m_FilFza = dofElemSet';
        m_Fil = reshape(repmat(reshape(dofElemSet,dofpe,[]),dofpe,1),1,[]);
        m_Col = reshape(repmat(m_FilFza,dofpe,1),1,[]);
        %Vector de incognitas al tiempo  "n+1" (Desplazamientos y poro
        %presiones)
        uElemSet  = reshape(u(dofElemSet),[],nElem);
        %Vector de incognitas al tiempo  "n" (Desplazamientos y poro
        %presiones)
        u_oldElemSet = reshape(u_old(dofElemSet),[],nElem);
    end
    %############################################################################################
    
    %Solo para debug, se guarda el numero de set en e_VG (ya que esta variable se pasa a todas las
    %funciones).
    e_VG.iSet = iSet;

    %Se podria poner el parfor de dentro de cada eltype, asi evitar hacer la verificacion de tipo
    %de elemento para cada set.
    %parfor iElem = 1:nElem
    %for iElem = 1:nElem

    %Esta linea no puede estar si esta el parfor activado.
    %e_VG.iElem = iElem;
    %No se puede modificar una variable definida previa al loop, si no es sliced y si no es
    %interpretada como de reduccion, ya que no sabe como interpretarla el MatLab (puede haber
    %superposicion de resultados, al hacer reduccion). Por eso cada Lab debe tener su "copia",
    %para modificarla en forma independiente.
    %Lo que puede ser lento es copiar para cada lab esa e_VG_Aux, que puede ser grande.
    %En realidad como cada procesador tiene su copia local del e_VG (ya que MatLab realiza una
    %copia por cada Lab de todas la variables, excepto las sliced, donde solo copia la parte
    %que le corresponde al procesador), y al "copiarse" esta al e_VG_Aux, lo que unico que se
    %hace es copiarse el puntero, ya que no se esta modificando el e_VG_Aux. Luego lo que se
    %realiza es la modificacion de un valor de un campo (field), donde supuestamente
    %MatLab al cambiar un field de una estructura no hace copia de toda la estructura de nuevo,
    %si solo del field, por lo tanto no tendria que ser mucho mas lenta.
    %Otra seria pasar la variable iElem como argumento en las funciones que llama dentro del
    %parfor.    
    
    
%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
%Bucle por elemento. P/calculo de fuerza interna y tensor tangente del elemento
%     parfor iElem = 1:nElem
    for iElem = 1:nElem
        %Estructura de variables globales
        e_VG_Aux16 = e_VG;
        %Numero de elemento
        e_VG_Aux16.iElemSet = iElem;
        e_VG_Aux16.iElemNum = m_NumElem(iElem);
        %Matriz de coordenadas
        coord_n = f_CoordElem(xx,conec(iElem,:));
        %Cuadrangulo estandar de 8 nodos (Aplicado p/medio bifase)    
        %Funcion que permite determinar matriz tangente, fuerza interna,
        %variables primales y duales y operadores tangenres
        [m_Ke(:,:,iElem),m_Fint(:,iElem),eps_new(:,iElem),phi_new(:,iElem),...
            porpr_new(:,iElem),eps_fluct(:,iElem),phi_fluct(:,iElem),p_fluct(:,iElem),...
            sigmaE_new(:,iElem),sigmaT_new(:,iElem),mflu_new(:,iElem),...
            velflu_sta(:,iElem),velflu_total(:,iElem),velflu_new(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
            m_sig_eps(:,:,:,iElem),m_chi_eps(:,:,:,iElem),m_V_eps(:,:,:,iElem),...
            m_sig_p(:,:,:,iElem),m_chi_p(:,:,:,iElem),m_V_p(:,:,:,iElem),....
            m_sig_phi(:,:,:,iElem),m_chi_phi(:,:,:,iElem),m_V_phi(:,:,:,iElem)]=...
            f_MatElem_BifaseMulTSc(uElemSet(:,iElem),u_oldElemSet(:,iElem),...
            coord_n,hvar_old(:,iElem),aux_var(:,iElem),m_BT_d(:,:,:,iElem),...
            m_DetJT_d(:,iElem),m_DerCa_p(:,:,:,iElem),m_DetJT_p(:,iElem),...
            m_FF_d,m_FF_p,m_DefMacroSet(:,iElem),m_GradPorMacroSet_dup(:,iElem),...
            m_GradPorMacroSet_up_n(:,iElem),m_PorMacroSet(:,iElem),m_DefMacroSet_new(:,iElem),...
            m_GradPorMacroSet_up_new(:,iElem),m_PorMacroSet_new(:,iElem),e_DatElemSet,e_DatMatSet,e_VG_Aux16);
        
    end %for(iElem) 
      
    %############################################################################################
    %SALIDA DE INFORMACION
    if e_VG.conshyp==50 || e_VG.conshyp==60 || e_VG.conshyp==61 || e_VG.conshyp==62  || (~isfield(e_VG,'conshypMacro') && e_VG.conshyp==14)
        %############################################################################################
        %VARIABLES MACRO-ESCALA
        %############################################################################################
        %VARIABLES MACRO-ESCALA EN EL PASO DE TIEMPO "n+1"
        e_VarEst_new(iSet) = f_OutVar_time(e_VarEst_new(iSet),e_VarEst_old(iSet),...
            eps_new,phi_new,porpr_new,eps_fluct,phi_fluct,p_fluct,sigmaE_new,sigmaT_new,...
            mflu_new,velflu_sta,velflu_total,velflu_new,e_VG,0);
        %############################################################################################
        
    elseif isfield(e_VG,'conshypMacro') 
        %PARA EL PASO MICRO - MACRO
        %############################################################################################
        %VARIABLES MICRO-ESCALA 
        e_VarEst_new(iSet) = f_OutVar_time(e_VarEst_new(iSet),e_VarEst_old(iSet),...
            eps_new,phi_new,porpr_new,eps_fluct,phi_fluct,p_fluct,sigmaE_new,sigmaT_new,...
            mflu_new,velflu_sta,velflu_total,velflu_new,e_VG,2);
        %############################################################################################        
    end
    %Variables historicas de cada modelo constitutivo al tiempo "n+1"   
    e_VarEst_new(iSet).hvar = hvar_new;
    e_VarEst_new(iSet).VarHistElem = m_VarHistElemNew;
    e_VarAux(iSet).VarAuxGP = aux_var;
    e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
    
    %Matriz tangente
    c_Ke{iSet} = m_Ke(:);
    %Vector de fuerzas internas
    c_Fint{iSet} = m_Fint(:);
    
    c_Fil{iSet} = m_Fil;
    c_Col{iSet} = m_Col;
    c_FilFza{iSet} = m_FilFza;

    %OPERADORES TANGENTES
    %Operador tangente constitutivo dsigma/deps
    c_sig_eps{iSet} = m_sig_eps;
    %Operador tangente constitutivo dchi/deps
    c_chi_eps{iSet}= m_chi_eps;
    %Operador tangente constitutivo dV/deps
    c_V_eps{iSet}= m_V_eps;
    %Operador tangenteconstitutivo dsigma/dp
    c_sig_p{iSet}= m_sig_p;
    %Operador tangente constitutivo dchi/dp    
    c_chi_p{iSet}= m_chi_p;
    %Operador tangente constitutivo dV/dp       
    c_V_p{iSet}= m_V_p;
    %Operador tangente constitutivo dsigma/dphi    
    c_sig_phi{iSet}= m_sig_phi;
    %Operador tangente constitutivo dchi/dphi    
    c_chi_phi{iSet}= m_chi_phi;
    %Operador tangente constitutivo dV/dphi    
    c_V_phi{iSet}= m_V_phi;
    
end %for(iSet)

% Ensamble de matriz de fuerzas internas global
Fint = sparse([c_FilFza{:}],1,cat(1,c_Fint{:}));

% Ensamble de matriz de rigidez tangente global
KT = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:}));

end
