function    [e_TanOp,m_VarFluc_eps0,m_VarFluc_p0,m_VarFluc_phi0] ...
    = f_OpTangHomBifase(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,Omega_micro,...
   m_ElemLocHomog,m_ElemLocPert,xx,e_VG,m_VarFluc_eps0,m_VarFluc_p0,m_VarFluc_phi0)

   %Ver que si no conviene por velocidad tener un funcion separada cuando se homogeniza en el
   %dominio localizado.
   
   % DETERMINACION DE LOS OPERADORES TANGENTES HOMOGENEIZADOS
   ntens = e_VG.ntens;
   nSet = e_VG.nSet;
   Dtime = e_VG.Dtime; 
   theta = e_VG.theta; 
  
   % Variacion de las fluctuaciones de los desplazamientos y poropresiones micros
   %Creo varias operaciones pueden ser evitadas, si se junta la siguiente funcion con esta, 
   %principalmente si se vectorizan las operaciones (ambas son integraciones sobre los elementos, 
   %pero que se debe realizar en forma separada, ya que hay resolver un sistema para obtener primero
   %los incrementos de las fluctuaciones m_VarFluc).
   [m_VarFluc_eps] = f_VarFluct_eps(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,m_ElemLocPert,e_VG,m_VarFluc_eps0);
   [m_VarFluc_p] = f_VarFluct_p(KT,m_LinCond,dofl,doff,e_DatSet,m_ElemLocPert,e_VG,m_VarFluc_p0);
   [m_VarFluc_phi] = f_VarFluct_phi(KT,m_LinCond,dofl,doff,e_DatSet,xx,m_ElemLocPert,e_VG,m_VarFluc_phi0);
   m_VarFluc_eps0 = m_VarFluc_eps0+m_VarFluc_eps;
   m_VarFluc_p0 = m_VarFluc_p0+m_VarFluc_p;
   m_VarFluc_phi0 = m_VarFluc_phi0+m_VarFluc_phi;
   % Homogenizacion
   % Operadores tangentes respecto de las deformaciones (eps)
   m_Csig_eps = zeros(ntens,ntens);
   m_bw_eps = zeros(1,ntens);
   m_Cw_eps = zeros(2,ntens);
   m_Cw_eps2 = zeros(2,ntens);
   m_Cw_eps3 = zeros(2,ntens);
   % Operadores tangentes respecto de las poropresiones (p)
   m_bsig_p = zeros(ntens,1);
   m_Mww_p = zeros(1,1);
   m_bw_p = zeros(2,1);
   m_bw_p2 = zeros(2,1);
   m_bw_p3 = zeros(2,1);
   % Operadores tangentes respecto del gradiente de las poropresiones (phi)
   m_ksig_phi = zeros(ntens,2);
   m_kww_phi = zeros(1,2);
   m_kw_phi = zeros(2,2);
   m_kw_phi2 = zeros(2,2);
   m_kw_phi3 = zeros(2,2);
   %Matriz identidad tensorial (se usa para ver si acelera el calculo, al sacar como factor comun
   %m_CT(:,:,iPG,iElem)). AA: D(hom)=D(Taylor)+D(fluctuante). Ambos son
   %funciones de D(micro)= m_CT entonces saca factor comun al calcular m_CTHomog
   m_IndTens = eye(ntens);
for iSet = 1:nSet
    
    nElem = e_DatSet(iSet).nElem;
    %conec = e_DatSet(iSet).conec;
    m_DofElem = e_DatSet(iSet).m_DofElem;
    m_DetJT = e_DatSet(iSet).m_DetJT_d;
    m_BT = e_DatSet(iSet).m_BT_d;
    m_CT = c_CT{iSet};
    e_DatElemSet = e_DatSet(iSet).e_DatElem;
    npg = e_DatElemSet.npg;
    wg = e_DatElemSet.wg;
    N4 = e_DatSet(iSet).m_FF_p;
    m_B_p = e_DatSet(iSet).m_DerCa_p;
    pos_d =  e_DatElemSet.pos_d;
    pos_p =  e_DatElemSet.pos_p;
    BiotM = e_DatSet(iSet).e_DatMat.m_Biot;
    beta = e_DatSet(iSet).e_DatMat.beta;
    PermK = e_DatSet(iSet).e_DatMat.m_PermK;
    conec = e_DatSet(iSet).conec;
    m_ElemLocSet = m_ElemLocHomog(e_DatSet(iSet).m_IndElemSet);
    for iElem = 1:nElem
        yy = f_CoordElem(xx,conec(iElem,:));
        %Se ensambla el elemento si esta en el dominio localizado.
        if m_ElemLocSet(iElem)
            %dofElem = f_DofElem(conec(iElem,:),ndn);
            dofElem = m_DofElem(:,iElem);
            m_pesoPG = m_DetJT(:,iElem).*wg;
            
            m_VarFlucElem_eps_u = m_VarFluc_eps(dofElem(pos_d),:);
            m_VarFlucElem_eps_p =  m_VarFluc_eps(dofElem(pos_p),:);
            
            m_VarFlucElem_eps0_p =  m_VarFluc_eps0(dofElem(pos_p),:);
            
            m_VarFlucElem_p_u = m_VarFluc_p(dofElem(pos_d),:);
            m_VarFlucElem_p_p = m_VarFluc_p(dofElem(pos_p),:);
            
            m_VarFlucElem_p0_p = m_VarFluc_p0(dofElem(pos_p),:);
            
            m_VarFlucElem_phi_u= m_VarFluc_phi(dofElem(pos_d),:);
            m_VarFlucElem_phi_p= m_VarFluc_phi(dofElem(pos_p),:);
            
            m_VarFlucElem_phi0_p= m_VarFluc_phi0(dofElem(pos_p),:);
            
            for iPG = 1:npg
                % %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % Operadores tangentes respecto de las deformaciones (eps)
                m_Csig_eps = m_Csig_eps+(m_CT(:,:,iPG,iElem)*(m_IndTens+...
                m_BT(:,:,iPG,iElem)*m_VarFlucElem_eps_u)-...
                BiotM*N4(:,:,iPG)*m_VarFlucElem_eps_p)*m_pesoPG(iPG);

                m_bw_eps = m_bw_eps+(BiotM'*(m_IndTens+...
                m_BT(:,:,iPG,iElem)*m_VarFlucElem_eps_u)+...
                beta*N4(:,:,iPG)*m_VarFlucElem_eps_p)*m_pesoPG(iPG);

                m_Cw_eps3 = m_Cw_eps3+(-PermK*m_B_p(:,:,iPG,iElem)*m_VarFlucElem_eps_p*m_pesoPG(iPG));

                m_Cw_eps2 = m_Cw_eps2-(((BiotM'*(m_IndTens+m_BT(:,:,iPG,iElem)*m_VarFlucElem_eps_u))'*...
                N4(:,:,iPG)*(yy(1:2,1:4))')'+...
                ((beta*N4(:,:,iPG)*m_VarFlucElem_eps_p)'*N4(:,:,iPG)*(yy(1:2,1:4))')')*m_pesoPG(iPG);  

                m_Cw_eps = m_Cw_eps+(-PermK*m_B_p(:,:,iPG,iElem)*m_VarFlucElem_eps0_p*m_pesoPG(iPG));

                % Operadores tangentes respecto de las poropresiones (p)
                m_bsig_p = m_bsig_p+(m_CT(:,:,iPG,iElem)*m_BT(:,:,iPG,iElem)*m_VarFlucElem_p_u-...
                BiotM*(1+N4(:,:,iPG)*m_VarFlucElem_p_p))*m_pesoPG(iPG);

                m_Mww_p =  m_Mww_p+((BiotM'*m_BT(:,:,iPG,iElem)*m_VarFlucElem_p_u)+...
                beta*(1+N4(:,:,iPG)*m_VarFlucElem_p_p))*m_pesoPG(iPG);

                m_bw_p3 = m_bw_p3+(-PermK*m_B_p(:,:,iPG,iElem)*m_VarFlucElem_p_p*m_pesoPG(iPG));

                m_bw_p2 = m_bw_p2-((BiotM'*m_BT(:,:,iPG,iElem)*m_VarFlucElem_p_u*N4(:,:,iPG)*(yy(1:2,1:4))'+...
                beta*(1+N4(:,:,iPG)*m_VarFlucElem_p_p)*N4(:,:,iPG)*(yy(1:2,1:4))')')*m_pesoPG(iPG);

                m_bw_p = m_bw_p+(-PermK*m_B_p(:,:,iPG,iElem)*m_VarFlucElem_p0_p*m_pesoPG(iPG));
                % Operadores tangentes respecto del gradiente de las poropresiones (phi)
                m_ksig_phi = m_ksig_phi+(m_CT(:,:,iPG,iElem)*m_BT(:,:,iPG,iElem)*m_VarFlucElem_phi_u-...
                BiotM*(N4(:,:,iPG)*(yy(1:2,1:4))'+N4(:,:,iPG)*m_VarFlucElem_phi_p))*m_pesoPG(iPG);

                m_kww_phi = m_kww_phi+(BiotM'*m_BT(:,:,iPG,iElem)*m_VarFlucElem_phi_u+...
                beta*(N4(:,:,iPG)*(yy(1:2,1:4))'+N4(:,:,iPG)*m_VarFlucElem_phi_p))*m_pesoPG(iPG);               

                m_kw_phi3 = m_kw_phi3+(-PermK*(eye(2)+m_B_p(:,:,iPG,iElem)*m_VarFlucElem_phi_p)*m_pesoPG(iPG)); 
                m_kw_phi2 = m_kw_phi2-((BiotM'*m_BT(:,:,iPG,iElem)*m_VarFlucElem_phi_u)'*N4(:,:,iPG)*(yy(1:2,1:4))'+...
                (beta*(N4(:,:,iPG)*(yy(1:2,1:4))'+...
                N4(:,:,iPG)*m_VarFlucElem_phi_p))'*N4(:,:,iPG)*(yy(1:2,1:4)'))*m_pesoPG(iPG); 

                m_kw_phi = m_kw_phi+(-PermK*(eye(2)+m_B_p(:,:,iPG,iElem)*m_VarFlucElem_phi0_p)*m_pesoPG(iPG)); 
                %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                          
            end
         end
    end      
end
   
   m_Csig_eps = m_Csig_eps/Omega_micro;
   m_bw_eps = m_bw_eps/Omega_micro;
   m_Cw_eps3 = m_Cw_eps3*theta*Dtime/Omega_micro;
   m_Cw_eps2 = m_Cw_eps2/Omega_micro;
   m_Cw_eps3 = m_Cw_eps3+m_Cw_eps2;
   
   m_Cw_eps = m_Cw_eps/Omega_micro;
   
   
   m_bsig_p = m_bsig_p/Omega_micro;
   m_Mww_p = m_Mww_p/Omega_micro;
   m_bw_p3 = m_bw_p3*theta*Dtime/Omega_micro;
   m_bw_p2 = m_bw_p2/Omega_micro;
   m_bw_p3 = m_bw_p3+m_bw_p2;
   
   m_bw_p = m_bw_p/Omega_micro;
      
   m_ksig_phi = m_ksig_phi/Omega_micro;
   m_kww_phi = m_kww_phi/Omega_micro;
   m_kw_phi3 = m_kw_phi3*theta*Dtime/Omega_micro;
   m_kw_phi2 = m_kw_phi2/Omega_micro;
   m_kw_phi3 = m_kw_phi3+m_kw_phi2;
  
   m_kw_phi = m_kw_phi/Omega_micro;

   e_TanOp = struct('m_Csig_eps',m_Csig_eps,'m_bw_eps',m_bw_eps,'m_Cw_eps3',m_Cw_eps3,...
                    'm_bsig_p',m_bsig_p,'m_Mww_p',m_Mww_p,'m_bw_p3',m_bw_p3,...
                    'm_ksig_phi',m_ksig_phi,'m_kww_phi',m_kww_phi,'m_kw_phi3',m_kw_phi3,...
                    'm_Cw_eps',m_Cw_eps,'m_bw_p',m_bw_p,'m_kw_phi',m_kw_phi);
end

%%
function [m_VarFluc]= ...
    f_VarFluct_eps(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,m_ElemLocPert,e_VG,m_VarFluc_eps0)
   %Numero de grados de libertad por nodo.
   %ndn = e_VG.ndn;
   ndoft = e_VG.ndoft;
   ntens = e_VG.ntens;
   nSet = e_VG.nSet;
   conshyp = e_VG.conshyp;
   % Representa derivadas respecto de tensor de 2do orden
   m_FCHom = zeros(ndoft,ntens); 
   Fext = zeros(ndoft,ntens); 

   for iSet = 1:nSet
      nElem = e_DatSet(iSet).nElem;
      m_DofElem = e_DatSet(iSet).m_DofElem;
      m_DetJT_d = e_DatSet(iSet).m_DetJT_d;
      m_BT = e_DatSet(iSet).m_BT_d;
      m_CT = c_CT{iSet};
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      npg = e_DatElemSet.npg;
      wg = e_DatElemSet.wg;
      dofpe_d = e_DatElemSet.dofpe_d;
      pos_d =  e_DatElemSet.pos_d;
      
      dofpe_p = e_DatElemSet.dofpe_p;
      N4 = e_DatSet(iSet).m_FF_p;
      m_DetJT_p = e_DatSet(iSet).m_DetJT_p;
      pos_p =  e_DatElemSet.pos_p;
      
      BiotM = e_DatSet(iSet).e_DatMat.m_Biot;


      m_ElemLocSet = m_ElemLocPert(e_DatSet(iSet).m_IndElemSet);

      for iElem = 1:nElem        
         
         %Se ensambla el elemento si esta en el dominio localizado.
         if m_ElemLocSet(iElem)
            m_pesoPG_d = m_DetJT_d(:,iElem).*wg;
            m_pesoPG_p = m_DetJT_p(:,iElem).*wg;

            m_FCHomElem_Ct = zeros(dofpe_d,ntens);
            m_FCHomElem_Biot = zeros(dofpe_p,ntens);
            for iPG = 1:npg
               %No se coloca el signo menos porque en incremental_disp cuando se resuelve el sistema ya
               %se lo pone.
               m_FCHomElem_Ct = m_FCHomElem_Ct+m_BT(:,:,iPG,iElem)'*m_CT(:,:,iPG,iElem)*m_pesoPG_d(iPG);
               %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               m_FCHomElem_Biot = m_FCHomElem_Biot-N4(:,:,iPG)'*BiotM'*m_pesoPG_p(iPG); %POR INTEGRACION TEMPORAL
               %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                m_FCHomElem_Biot = m_FCHomElem_Biot+N4(:,:,iPG)'*BiotM'*m_pesoPG_p(iPG);
%                m_FCHomElem_Biot = m_FCHomElem_Biot+N4(:,:,iPG,iElem)'*BiotM'*m_pesoPG_p(iPG);
            end
            %dofElem = f_DofElem(conec(iElem,:),ndn);
            dofElem = m_DofElem(:,iElem);
            m_FCHom(dofElem(pos_d),:) = m_FCHom(dofElem(pos_d),:)+m_FCHomElem_Ct;
            m_FCHom(dofElem(pos_p),:) = m_FCHom(dofElem(pos_p),:)+m_FCHomElem_Biot;
         end
         
      end
      
   end
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if conshyp == 15
       [Fext] = f_Fflux2(Fext,m_VarFluc_eps0,e_DatSet,e_VG); %AA: add function
   elseif conshyp ==16
       [Fext] = f_Fflux_ML2(Fext,m_VarFluc_eps0,e_DatSet,e_VG); %AA: add function
   end
   m_FCHom = m_FCHom - Fext;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   %No se usa la funcion residuo porque no maneja el caso del vector de fuerzas con muchas columnas.
   %Este problema es directo, es decir no se esta calculando diferencia de fuerzas (residuo), sino
   %que son "fuerzas" aplicadas al problema (condiciones de borde de tipo fuerzas).
   m_FCHom = m_FCHom(dofl,:)+m_LinCond'*m_FCHom(doff,:);
   %Determinacion del incremento o variacion de las fluctuaciones.
   m_VarFluc = zeros(ndoft,ntens);
   m_VarFluc(dofl,:) = incremental_disp(m_FCHom,KT,0,m_LinCond,dofl,doff,[],[],e_VG);
   %En los grados de libertad restringidos no hay variacion de las fluctuaciones, es decir como si
   %las condiciones de borde constantes para este problema fuera siempre nulas. No aporta a ecuacion
   %siguiente al ser C = 0 (uR = L*uP+C).
   m_VarFluc(doff,:) = m_LinCond*m_VarFluc(dofl,:);   
   
end
%%
function [m_VarFluc] = ...
    f_VarFluct_p(KT,m_LinCond,dofl,doff,e_DatSet,m_ElemLocPert,e_VG,m_VF_p_old)

   ndoft = e_VG.ndoft;
   nSet = e_VG.nSet;
   conshyp = e_VG.conshyp;
   
   m_FCHom = zeros(ndoft,1); % Representa derivadas respecto de un escalar
   Fext = zeros(ndoft,1); 
   
   for iSet = 1:nSet
      nElem = e_DatSet(iSet).nElem;
      m_DofElem = e_DatSet(iSet).m_DofElem;
      m_DetJT_d = e_DatSet(iSet).m_DetJT_d;
      m_BT = e_DatSet(iSet).m_BT_d;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      npg = e_DatElemSet.npg;
      wg = e_DatElemSet.wg;
      dofpe_d = e_DatElemSet.dofpe_d;
      pos_d =  e_DatElemSet.pos_d;
      
      dofpe_p = e_DatElemSet.dofpe_p;
      N4 = e_DatSet(iSet).m_FF_p;
      m_DetJT_p = e_DatSet(iSet).m_DetJT_p;
      pos_p =  e_DatElemSet.pos_p;
      
      BiotM = e_DatSet(iSet).e_DatMat.m_Biot;
      beta = e_DatSet(iSet).e_DatMat.beta;


      m_ElemLocSet = m_ElemLocPert(e_DatSet(iSet).m_IndElemSet);

      for iElem = 1:nElem        
         
         %Se ensambla el elemento si esta en el dominio localizado.
         if m_ElemLocSet(iElem)
            m_pesoPG_d = m_DetJT_d(:,iElem).*wg;
            m_pesoPG_p = m_DetJT_p(:,iElem).*wg;

            m_FCHomElem_Bi = zeros(dofpe_d,1);
            m_FCHomElem_BiotMod = zeros(dofpe_p,1);
            for iPG = 1:npg
               %No se coloca el signo menos porque en incremental_disp cuando se resuelve el sistema ya
               %se lo pone.
               m_FCHomElem_Bi = m_FCHomElem_Bi-m_BT(:,:,iPG,iElem)'*BiotM*m_pesoPG_d(iPG);
               %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               m_FCHomElem_BiotMod = m_FCHomElem_BiotMod-N4(:,:,iPG)'*beta*m_pesoPG_p(iPG);
               %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                m_FCHomElem_BiotMod = m_FCHomElem_BiotMod+N4(:,:,iPG,iElem)'*beta*m_pesoPG_p(iPG);
            end
            %dofElem = f_DofElem(conec(iElem,:),ndn);
            dofElem = m_DofElem(:,iElem);
            m_FCHom(dofElem(pos_d),:) = m_FCHom(dofElem(pos_d),:)+m_FCHomElem_Bi;
            m_FCHom(dofElem(pos_p),:) = m_FCHom(dofElem(pos_p),:)+m_FCHomElem_BiotMod;
         end
         
      end
      
   end
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if conshyp == 15
       [Fext] = f_Fflux(Fext,m_VF_p_old,e_DatSet,e_VG); %AA: add function
   elseif conshyp == 15
       [Fext] = f_Fflux_ML(Fext,m_VF_p_old,e_DatSet,e_VG); %AA: add function
   end
   m_FCHom = m_FCHom - Fext;
  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   %No se usa la funcion residuo porque no maneja el caso del vector de fuerzas con muchas columnas.
   %Este problema es directo, es decir no se esta calculando diferencia de fuerzas (residuo), sino
   %que son "fuerzas" aplicadas al problema (condiciones de borde de tipo fuerzas).
   m_FCHom = m_FCHom(dofl,:)+m_LinCond'*m_FCHom(doff,:);  
   %Determinacion del incremento o variacion de las fluctuaciones.
   m_VarFluc = zeros(ndoft,1);
   m_VarFluc(dofl,:) = incremental_disp(m_FCHom,KT,0,m_LinCond,dofl,doff,[],[],e_VG);
   %En los grados de libertad restringidos no hay variacion de las fluctuaciones, es decir como si
   %las condiciones de borde constantes para este problema fuera siempre nulas. No aporta a ecuacion
   %siguiente al ser C = 0 (uR = L*uP+C).
   m_VarFluc(doff,:) = m_LinCond*m_VarFluc(dofl,:);   
   
end
%%
function [m_VarFluc] = ...
    f_VarFluct_phi(KT,m_LinCond,dofl,doff,e_DatSet,xx,m_ElemLocPert,e_VG,m_VF_phi_old)

   ndoft = e_VG.ndoft;
   nSet = e_VG.nSet;
   conshyp = e_VG.conshyp;
   
   m_FCHom = zeros(ndoft,2); % Representa derivadas respecto de un vector
   m_FCHom_theta = zeros(ndoft,2); % Representa derivadas respecto de un vector
   Fext = zeros(ndoft,2); 
   for iSet = 1:nSet
      nElem = e_DatSet(iSet).nElem;
      conec = e_DatSet(iSet).conec;
      m_DofElem = e_DatSet(iSet).m_DofElem;
      m_DetJT_d = e_DatSet(iSet).m_DetJT_d;
      m_BT = e_DatSet(iSet).m_BT_d;
      N4 = e_DatSet(iSet).m_FF_p;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      npg = e_DatElemSet.npg;
      wg = e_DatElemSet.wg;
      dofpe_d = e_DatElemSet.dofpe_d;
      pos_d =  e_DatElemSet.pos_d;
      
      dofpe_p = e_DatElemSet.dofpe_p;
      m_B_p = e_DatSet(iSet).m_DerCa_p;
      m_DetJT_p = e_DatSet(iSet).m_DetJT_p;
      pos_p =  e_DatElemSet.pos_p;
      
      BiotM = e_DatSet(iSet).e_DatMat.m_Biot;
      beta = e_DatSet(iSet).e_DatMat.beta;
      PermK = e_DatSet(iSet).e_DatMat.m_PermK;

      m_ElemLocSet = m_ElemLocPert(e_DatSet(iSet).m_IndElemSet);

      for iElem = 1:nElem        
         yy = f_CoordElem(xx,conec(iElem,:));
         %Se ensambla el elemento si esta en el dominio localizado.
         if m_ElemLocSet(iElem)
            m_pesoPG_d = m_DetJT_d(:,iElem).*wg;
            m_pesoPG_p = m_DetJT_p(:,iElem).*wg;

            m_FCHomElem_Biy = zeros(dofpe_d,2);
            m_FCHomElem_Per1 = zeros(dofpe_p,2);
            m_FCHomElem_Per2 = zeros(dofpe_p,2);
            for iPG = 1:npg
               %No se coloca el signo menos porque en incremental_disp cuando se resuelve el sistema ya
               %se lo pone.
               m_FCHomElem_Biy = m_FCHomElem_Biy-(m_BT(:,:,iPG,iElem)'*BiotM*...
                   N4(:,:,iPG)*(yy(1:2,1:4))')*m_pesoPG_d(iPG);
               
               m_FCHomElem_Per1 = m_FCHomElem_Per1-(beta*N4(:,:,iPG)'*N4(:,:,iPG)*(yy(1:2,1:4))')*m_pesoPG_p(iPG);
               
               m_FCHomElem_Per2 = m_FCHomElem_Per2-(m_B_p(:,:,iPG,iElem)'*PermK)*m_pesoPG_p(iPG);
            end
            %dofElem = f_DofElem(conec(iElem,:),ndn);
            dofElem = m_DofElem(:,iElem);
            m_FCHom(dofElem(pos_d),:) = m_FCHom(dofElem(pos_d),:)+m_FCHomElem_Biy;
            m_FCHom(dofElem(pos_p),:) = m_FCHom(dofElem(pos_p),:)+m_FCHomElem_Per1;
            m_FCHom_theta(dofElem(pos_p),:) = m_FCHom_theta(dofElem(pos_p),:)+m_FCHomElem_Per2;
         end
         
      end
      
   end
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   m_FCHom(3:3:end,:) = m_FCHom(3:3:end,:) + m_FCHom_theta(3:3:end,:)*e_VG.Dtime*(e_VG.theta+e_VG.istep);
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if conshyp == 15
       [Fext] = f_Fflux3(Fext,m_VF_phi_old,e_DatSet,e_VG); %AA: add function
   elseif conshyp ==16
       [Fext] = f_Fflux_ML3(Fext,m_VF_phi_old,e_DatSet,e_VG); %AA: add function
   end
   m_FCHom = m_FCHom - Fext;
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%    m_FCHom = sparse(m_Fil(:),m_Col(:),m_FCHom(:),ndoft,ntens);   
   %No se usa la funcion residuo porque no maneja el caso del vector de fuerzas con muchas columnas.
   %Este problema es directo, es decir no se esta calculando diferencia de fuerzas (residuo), sino
   %que son "fuerzas" aplicadas al problema (condiciones de borde de tipo fuerzas).
   %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   m_FCHom = m_FCHom(dofl,:)+m_LinCond'*m_FCHom(doff,:);
   %Determinacion del incremento o variacion de las fluctuaciones.
   m_VarFluc = zeros(ndoft,2);
   m_VarFluc(dofl,:) = incremental_disp(m_FCHom,KT,0,m_LinCond,dofl,doff,[],[],e_VG);
   %En los grados de libertad restringidos no hay variacion de las fluctuaciones, es decir como si
   %las condiciones de borde constantes para este problema fuera siempre nulas. No aporta a ecuacion
   %siguiente al ser C = 0 (uR = L*uP+C).
   m_VarFluc(doff,:) = m_LinCond*m_VarFluc(dofl,:);   
   
end
