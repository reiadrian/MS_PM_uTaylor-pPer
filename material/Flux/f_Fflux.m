function [Fext] = f_Fflux(Fext,u,e_DatSet,e_VG)
% Determina las fuerzas debidado a las poropresiones en cada paso de tiempo
% y los suma al vector global de fuerzas externas
%Sirve tanto para caso mono-escala como para la micro-escala
%de un caso multiescala

nSet = e_VG.nSet;
c_Fext = cell(nSet,1);
c_FilFza = cell(nSet,1);

for iSet = 1:nSet
    
    Dtime = e_VG.Dtime;
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
    uElemSet  = reshape(u(dofElemSet),[],nElem);
    
    m_Fext = zeros(dofpe,nElem);
    for iElem = 1:nElem
        Hww = zeros(dofpe_p,dofpe_p);
        m_Dercae_p = m_DerCa_p(:,:,:,iElem);
        ue = uElemSet(:,iElem);
        ue_p = ue(pos_p);
        m_pesoPG_p = m_DetJT_p(:,iElem).*wg;
    for iPG = 1:nPG
        DerivN = m_Dercae_p(:,:,iPG);
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); 
    end %for(iPG)
    
    m_Fext(pos_p,iElem) =  Dtime*Hww*ue_p;
    m_Fext(pos_p0,iElem) =  0.0; %Es necesario?
    
    end %for(iElem)
        
    c_Fext{iSet} = m_Fext(:);
    c_FilFza{iSet} = m_FilFza;
end %for(iSet)
%###################################################################################
if e_VG.conshyp==17
    Fext(1:e_VG.ndoft-3) = Fext(1:e_VG.ndoft-3) + sparse([c_FilFza{:}],1,cat(1,c_Fext{:}));
else
    % Ensamble de matriz de fuerzas externas global
    Fext = Fext + sparse([c_FilFza{:}],1,cat(1,c_Fext{:}));
end
%###################################################################################


    
    

