function [Fext] = f_Fflux_phi_Pen(Fext,u,e_DatSet,e_VG)
%Determina las fuerzas debido a la relacion entre poro presiones fluctuantes
%y macro-gradientes de poro presiones en cada paso de tiempo
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
    m_DetJT_p = e_DatSet(iSet).m_DetJT_p;
    m_DofElem = e_DatSet(iSet).m_DofElem;
    nElem = e_DatSet(iSet).nElem;
     
    dofpe = e_DatElemSet.dofpe;
    dofpe_p = e_DatElemSet.dofpe_p;
    wg = e_DatElemSet.wg;
    nPG = e_DatElemSet.npg;
    pos_p = e_DatElemSet.pos_p;
    pos_p0 = e_DatElemSet.pos_p0;
    
    PermK = e_DatMatSet.m_PermK;
       
    % Grados de libertad y coordenadas de los nodos de los elementos del set
    dofElemSet = m_DofElem(:); 
    m_FilFza = dofElemSet';
    uElemSet_exx  = reshape(uxx(dofElemSet),[],nElem);
    uElemSet_eyy  = reshape(uyy(dofElemSet),[],nElem);
    m_Fextxx = zeros(dofpe,nElem);
    m_Fextyy = zeros(dofpe,nElem);
    
    for iElem = 1:nElem
        Hww = zeros(dofpe_p,dofpe_p);
        m_Dercae_p = m_DerCa_p(:,:,:,iElem);
        ue_exx = uElemSet_exx(:,iElem);
        ue_eyy = uElemSet_eyy(:,iElem);
        ue_p_exx = ue_exx(pos_p);
        ue_p_eyy = ue_eyy(pos_p);
        m_pesoPG_p = m_DetJT_p(:,iElem).*wg;
    for iPG = 1:nPG
        DerivN = m_Dercae_p(:,:,iPG);
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); 
    end %for(iPG)
    
    m_Fextxx(pos_p,iElem) =  Dtime*Hww*ue_p_exx;
    m_Fextxx(pos_p0,iElem) =  0.0; %Es necesario?
    
    m_Fextyy(pos_p,iElem) =  Dtime*Hww*ue_p_eyy;
    m_Fextyy(pos_p0,iElem) =  0.0; %Es necesario?
    
       
    end %for(iElem)
        
    c_Fextxx{iSet} = m_Fextxx(:);
    c_FilFzaxx{iSet} = m_FilFza;
    
    c_Fextyy{iSet} = m_Fextyy(:);
    c_FilFzayy{iSet} = m_FilFza;

end %for(iSet)
% Fext(1:3:end) = -Dtime*Fext(3:3:end); %ESTA BIEN?????????
% Ensamble de matriz de fuerzas externas global
%##########################################################################################
if e_VG.conshyp==17
    Fext(1:e_VG.ndoft-3,1) = Fext(1:e_VG.ndoft-3,1) + sparse([c_FilFzaxx{:}],1,cat(1,c_Fextxx{:}));
    Fext(1:e_VG.ndoft-3,2) = Fext(1:e_VG.ndoft-3,2) + sparse([c_FilFzayy{:}],1,cat(1,c_Fextyy{:}));
else
    % Ensamble de matriz de fuerzas externas global
    Fext(:,1) = Fext(:,1) + sparse([c_FilFzaxx{:}],1,cat(1,c_Fextxx{:}));
    Fext(:,2) = Fext(:,2) + sparse([c_FilFzayy{:}],1,cat(1,c_Fextyy{:}));
end
%###################################################################################
