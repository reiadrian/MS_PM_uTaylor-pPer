function [Fext] = f_Fflux_phi_ML(Fext,u,e_DatSet,e_VG)
%Determina las fuerzas debido a la relacion entre poro presiones fluctuantes
%y macro-gradientes de poro presiones
%en cada paso de tiempo considerando los terminos
%adicionales dados por los multiplicadores de Lagrange
%y los suma al vector global de fuerzas externas
%Matriz donde cada columna representa la fuerza ocasiona por:
%fx(p_fluct/phi_x)  fy(p_fluct/phi_y)

nSet = e_VG.nSet;
c_Fextxx = cell(nSet,1);
c_FilFzaxx = cell(nSet,1);
c_Fextyy = cell(nSet,1);
c_FilFzayy = cell(nSet,1);
Dtime = e_VG.Dtime;
uxx = u(:,1);
uyy = u(:,2);
for iSet = 1:nSet
    e_DatMatSet = e_DatSet(iSet).e_DatMat;
    e_DatElemSet = e_DatSet(iSet).e_DatElem;
    m_DerCa_p = e_DatSet(iSet).m_DerCa_p;
    m_FF_p = e_DatSet(iSet).m_FF_p;
    m_DetJT_p = e_DatSet(iSet).m_DetJT_p;
    m_DofElem = e_DatSet(iSet).m_DofElem;
    nElem = e_DatSet(iSet).nElem;
     
    dofpe = e_DatElemSet.dofpe;
    dofpe_p = e_DatElemSet.dofpe_p;
    wg = e_DatElemSet.wg;
    nPG = e_DatElemSet.npg;
    pos_p = e_DatElemSet.pos_p;
    pos_p0 = e_DatElemSet.pos_p0;
    pos_lambda = e_DatElemSet.pos_lambda;
%     pos_lambda0 = e_DatElemSet.pos_lambda0;
    
    PermK = e_DatMatSet.m_PermK;
       
    
       %###########################################################################
    if e_VG.conshyp==17
        ndoft=e_VG.ndoft;
                
        % Grados de libertad y coordenadas de los nodos de los elementos del set
        dofElemSet_p = [m_DofElem; repmat(ndoft,1,nElem)];
        dofElemSet = dofElemSet_p(:);
        m_FilFza = dofElemSet';
        
        
        uElemSet_exx  = reshape(uxx(dofElemSet),[],nElem);
        uElemSet_eyy  = reshape(uyy(dofElemSet),[],nElem);
    
    else
    
        % Grados de libertad y coordenadas de los nodos de los elementos del set
        dofElemSet = m_DofElem(:);
        m_FilFza = dofElemSet';
        
        uElemSet_exx  = reshape(uxx(dofElemSet),[],nElem);
        uElemSet_eyy  = reshape(uyy(dofElemSet),[],nElem);
    end
    %###########################################################################
    
    % Grados de libertad y coordenadas de los nodos de los elementos del set
%     dofElemSet = m_DofElem(:); 
%     m_FilFza = dofElemSet';
%     uElemSet_exx  = reshape(uxx(dofElemSet),[],nElem);
%     uElemSet_eyy  = reshape(uyy(dofElemSet),[],nElem);
    
    m_Fextxx = zeros(dofpe,nElem);
    m_Fextyy = zeros(dofpe,nElem);
    
    for iElem = 1:nElem
        Hww = zeros(dofpe_p,dofpe_p);
%         Lalfa = zeros(dofpe_p,dofpe_p);
        %Vector que contiene terminos adicionales por el Mult. Lagrange
        Lalfa = zeros(dofpe_p,1);
        m_Dercae_p = m_DerCa_p(:,:,:,iElem);
        ue_exx = uElemSet_exx(:,iElem);
        ue_eyy = uElemSet_eyy(:,iElem);
        ue_p_exx = ue_exx(pos_p);
        ue_p_eyy = ue_eyy(pos_p);
        
        ue_lambda_exx = ue_exx(pos_lambda);
        ue_lambda_eyy = ue_eyy(pos_lambda);
        
        m_pesoPG_p = m_DetJT_p(:,iElem).*wg;
    for iPG = 1:nPG
        N4 = m_FF_p(:,:,iPG);
        DerivN = m_Dercae_p(:,:,iPG);
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG);
%         Lalfa = Lalfa + N4'*N4*m_pesoPG_p(iPG);
        Lalfa = Lalfa + N4'*m_pesoPG_p(iPG); %Vector adicional para metodo de Mult. Lagrange
    end %for(iPG)
    
%     m_Fextxx(pos_p,iElem) =  Dtime*Hww*ue_p_exx;
    %ENTIENDO DEBERIA INCLUIR LA MATRIZ POR INCOGNITA ADICIONAL DE
    %MULTIPLOCADOR DE LAGRANGE PERO CUANDO LO HAGO EN DTIME ALTOS NO
    %CONVERGE
    %#######################################################################
    m_Fextxx(pos_p,iElem) =  Dtime*(Hww*ue_p_exx+Lalfa*ue_lambda_exx);
    m_Fextxx(pos_lambda,iElem) =  Dtime*Lalfa'*ue_p_exx;
    %#######################################################################
    m_Fextxx(pos_p0,iElem) =  0.0; %Es necesario?
%     m_Fextxx(pos_lambda0,iElem) =  0.0; %Es necesario?
    
%     m_Fextyy(pos_p,iElem) =  Dtime*Hww*ue_p_eyy;
    %ENTIENDO DEBERIA INCLUIR LA MATRIZ POR INCOGNITA ADICIONAL DE
    %MULTIPLOCADOR DE LAGRANGE PERO CUANDO LO HAGO EN DTIME ALTOS NO
    %CONVERGE
    %#######################################################################
    m_Fextyy(pos_p,iElem) =  Dtime*(Hww*ue_p_eyy+Lalfa*ue_lambda_eyy);
    m_Fextyy(pos_lambda,iElem) =  Dtime*Lalfa'*ue_p_eyy;
    %#######################################################################
    m_Fextyy(pos_p0,iElem) =  0.0; %Es necesario?
%     m_Fextyy(pos_lambda0,iElem) =  0.0; %Es necesario?
    
       
    end %for(iElem)
        
    c_Fextxx{iSet} = m_Fextxx(:);
    c_FilFzaxx{iSet} = m_FilFza;
    
    c_Fextyy{iSet} = m_Fextyy(:);
    c_FilFzayy{iSet} = m_FilFza;

end %for(iSet)
% Fext(1:3:end) = -Dtime*Fext(3:3:end); %ESTA BIEN?????????
% Ensamble de matriz de fuerzas externas global
Fext(:,1) = Fext(:,1) + sparse([c_FilFzaxx{:}],1,cat(1,c_Fextxx{:}));
Fext(:,2) = Fext(:,2) + sparse([c_FilFzayy{:}],1,cat(1,c_Fextyy{:}));