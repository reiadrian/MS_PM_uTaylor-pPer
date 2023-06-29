function Snapshots = SnapshotSave(Snapshots,istep,e_VarAux,e_VarEst_new,e_DatSet,e_VG,nglT)
% Funcion que guarda lo snapshot pedidos para cada paso equilibrado
% nota solo para: e_VG.isMICRO.MICRO=1
%                 eltype 108  %Cuadr�ngulo de 4 nodos FBar
% e_VG.Snap = vector binario que indica que snapshot debe guardarse

if ~e_VG.isMICRO.MICRO
    return
end

nset  = e_VG.nSet ;
ntens = e_VG.ntens ;
ind_flag = zeros(1,nset) ;

SnapStrain       = zeros(nglT*ntens,1) ;
SnapStress       = zeros(nglT*ntens,1) ;
SnapWeight       = zeros(nglT,1) ;
SnapEnergy_e     = zeros(nglT,1) ;
SnapEnergy_e_vol = zeros(nglT,1) ;
SnapEnergy_e_dev = zeros(nglT,1) ;
SnapEnergy_p     = zeros(nglT,1) ;
SnapEnergy_t     = zeros(nglT,1) ;
Snapflag         = 0 ;



for iset=1:nset    
    nElemS = e_DatSet(iset).nElem ; 
    npg    = e_DatSet(iset).e_DatElem.npg ; 
    
%     Psi_e     = zeros(npg,nElemS) ;
%     Psi_e_vol = zeros(npg,nElemS) ;
%     Psi_e_dev = zeros(npg,nElemS) ;
%     Psi_p     = zeros(npg,nElemS) ;
%     Psi_t     = zeros(npg,nElemS) ;

% Recupera valores de las variables historicas
    Strain_new = e_VarEst_new(iset).eps_fluct ;
    Stress_new = e_VarEst_new(iset).sigma ; % PK1
    
    ind  = 1+([1 2 3 4]-1)*npg + 1 ;
    
    Psi_e     = e_VarEst_new(iset).VarHistElem(ind,:) ;
    Psi_e_vol = e_VarEst_new(iset).VarHistElem(ind+1,:) ;
    Psi_e_dev = e_VarEst_new(iset).VarHistElem(ind+2,:) ;
    Psi_p     = e_VarEst_new(iset).VarHistElem(ind+3,:) ;
    Psi_t     = Psi_e + Psi_p ;
    
    Weight = e_DatSet(iset).m_DetJT ;% Falta multiplicar por el peso de gauss =1     

    Strain_new = reshape(Strain_new,[nElemS*npg*ntens 1]) ;
    Stress_new = reshape(Stress_new,[nElemS*npg*ntens 1]) ;
    Psi_e      = reshape(Psi_e     ,[nElemS*npg 1]) ;  
    Psi_e_vol  = reshape(Psi_e_vol ,[nElemS*npg 1]) ;
    Psi_e_dev  = reshape(Psi_e_dev ,[nElemS*npg 1]) ;
    Psi_p      = reshape(Psi_p     ,[nElemS*npg 1]) ;
    Psi_t      = reshape(Psi_t     ,[nElemS*npg 1]) ;
    Weight     = reshape(Weight    ,[nElemS*npg 1]) ;
        
    ElemSet = e_DatSet(iset).m_IndElemSet ; 
    ind_strain=[] ; ind_energy =[] ;        
    for ielem=1:nElemS
        ind_strain=[ind_strain ;[(ElemSet(ielem)-1)*npg*ntens + [1:1:npg*ntens]]'];  
        ind_energy=[ind_energy ;[(ElemSet(ielem)-1)*npg + [1:1:npg]]'];
    end
    SnapStrain(ind_strain,1) = Strain_new ; 
    SnapStress(ind_strain,1) = Stress_new ; 
    
%     if istep==1
        SnapEnergy_e(ind_energy,1)     = Psi_e ;
        SnapEnergy_e_vol(ind_energy,1) = Psi_e_vol ;
        SnapEnergy_e_dev(ind_energy,1) = Psi_e_dev ;
        SnapEnergy_p(ind_energy,1)     = Psi_p ;
        SnapEnergy_t(ind_energy,1)     = Psi_t ;
        SnapWeight(ind_energy,1)       = Weight ;
%     else
%         SnapEnergy_e(ind_energy,1)     = Snapshots.SnapEnergy_e(ind_energy,istep-1) + Psi_e ;
%         SnapEnergy_e_vol(ind_energy,1) = Snapshots.SnapEnergy_e_vol(ind_energy,istep-1) + Psi_e_vol ;
%         SnapEnergy_e_dev(ind_energy,1) = Snapshots.SnapEnergy_e_dev(ind_energy,istep-1) + Psi_e_dev ;
%         SnapEnergy_p(ind_energy,1)     = Snapshots.SnapEnergy_p(ind_energy,istep-1) + Psi_p ;
%         SnapEnergy_t(ind_energy,1)     = Snapshots.SnapEnergy_t(ind_energy,istep-1) + Psi_t ;
%     end 
    ind_flag(1,iset) = any(any(e_VarAux(iset).VarAuxGP,2),1) ;

end

Snapflag = any(ind_flag,2) ;

if e_VG.Snap(1)==1; Snapshots.SnapStrain(:,istep)       = SnapStrain ; end;
if e_VG.Snap(2)==1; Snapshots.SnapStress(:,istep)       = SnapStress ; end;
if istep==1 && e_VG.Snap(3)==1; Snapshots.SnapWeight    = SnapWeight ; end;
if istep==1
    if e_VG.Snap(4)==1; Snapshots.SnapEnergy_e(:,istep)     = SnapEnergy_e ; end;
    if e_VG.Snap(5)==1; Snapshots.SnapEnergy_e_vol(:,istep) = SnapEnergy_e_vol ; end;
    if e_VG.Snap(6)==1; Snapshots.SnapEnergy_e_dev(:,istep) = SnapEnergy_e_dev ; end;
    if e_VG.Snap(7)==1; Snapshots.SnapEnergy_p(:,istep)     = SnapEnergy_p ; end;
    if e_VG.Snap(8)==1; Snapshots.SnapEnergy_t(:,istep)     = SnapEnergy_t ; end;
    if e_VG.Snap(9)==1; Snapshots.Snapflag(1,istep)         = Snapflag ;  end;
else
    if e_VG.Snap(4)==1; 
        Snapshots.SnapEnergy_e(:,istep) = Snapshots.SnapEnergy_e(:,istep-1) + SnapEnergy_e ; 
    end;
    if e_VG.Snap(5)==1; 
        Snapshots.SnapEnergy_e_vol(:,istep) = Snapshots.SnapEnergy_e_vol(:,istep-1) + SnapEnergy_e_vol ; 
    end;
    if e_VG.Snap(6)==1; 
        Snapshots.SnapEnergy_e_dev(:,istep) = Snapshots.SnapEnergy_e_dev(:,istep-1) +  SnapEnergy_e_dev ; 
    end;
    if e_VG.Snap(7)==1; 
        Snapshots.SnapEnergy_p(:,istep) = Snapshots.SnapEnergy_p(:,istep-1) + SnapEnergy_p ; 
    end;
    if e_VG.Snap(8)==1; 
        Snapshots.SnapEnergy_t(:,istep) = Snapshots.SnapEnergy_t(:,istep-1) + SnapEnergy_t ; 
    end;
    if e_VG.Snap(9)==1; 
        Snapshots.Snapflag(:,istep) = Snapflag ; 
    end;
end


% save('Snapflag.mat'  ,'-struct', 'Snapshots', 'Snapflag');
% save('SnapStrain.mat','-struct', 'Snapshots', 'SnapStrain');
% save('SnapStress.mat','-struct', 'Snapshots', 'SnapStress');
% save('SnapWeight.mat','-struct', 'Snapshots', 'SnapWeight'); 
% % save('SnapEnergy.mat','Snapshots.SnapEnergy_e','Snapshots.SnapEnergy_e_vol',...
% %     'Snapshots.SnapEnergy_e_dev','Snapshots.SnapEnergy_p','Snapshots.SnapEnergy_t')
