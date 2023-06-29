function [kt,fint,sigmaE_new,sigmaT_new,eps_new,velflu_new,mflu_new,mfluY_new,...
    eps_fluct,phi_new,phi_fluct,pNodo_new,p_fluct,p_M,hvar_new,aux_var,m_Csig_eps,m_bw_eps,...
    m_Cw_eps3,m_bsig_p,m_Mww_p,m_bw_p3,m_ksig_phi,m_kww_phi,m_kw_phi3,m_Cw_eps,m_bw_p,m_kw_phi,m_TensorTang] = ...
   f_MatElem_BifaseMulTSc_old(u,eps_old,yy,hvar_old,aux_var,e_DatElemSet,e_DatMatSet,...
   m_Be_d,m_DetJe_d,m_Dercae_p,m_DetJe_p,m_FF_p,DefMacro,...
   GradPorMacro,PorMacro,sigmaE_old,sigmaT_old,e_VG)
%AA: Cree funcion
% Variable globales
ntens = e_VG.ntens;
theta = e_VG.theta;
Dtime = e_VG.Dtime;

dofpe = e_DatElemSet.dofpe;
dofpe_d = e_DatElemSet.dofpe_d; %AA
dofpe_p = e_DatElemSet.dofpe_p; %AA
nPG = e_DatElemSet.npg; 
wg = e_DatElemSet.wg;
pos_d = e_DatElemSet.pos_d;
pos_p = e_DatElemSet.pos_p;

sihvarpg = e_DatMatSet.sihvarpg;
siavarpg = e_DatMatSet.siavarpg;
conshyp = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;

% Inicializaciones
Kss = zeros(dofpe_d,dofpe_d);
Qsp = zeros(dofpe_d,dofpe_p);
Qps = zeros(dofpe_p,dofpe_d);
Sww = zeros(dofpe_p,dofpe_p);
Hww = zeros(dofpe_p,dofpe_p);

Kpp = zeros(dofpe_p,dofpe_p);

Lalfa = zeros(dofpe_p,dofpe_p);

% if isfield(e_VG,'conshypMacro') && conshyp == 15
%     alfa = e_DatMatSet.alfa_penalidad;  
% elseif isfield(e_VG,'conshypMacro') && conshyp == 16
%     fMul_Lag = zeros(dofpe_p,1); %Fuerzas internas asociadas al multiplicador de Lagrange
% end

fstreT = zeros(dofpe_d,1); %Fuerzas internas debidas a las tensiones totales

fflux_m = zeros(dofpe_p,1); %Fuerzas internas asociadas al flujo por el contenido de masa del fluido
fflux_V = zeros(dofpe_p,1); %Fuerzas internas asociadas al flujo por la velocidad de filtracion
fflux_mY = zeros(dofpe_p,1); %Fuerzas internas asociadas al flujo por la velocidad de filtracion
fflux = zeros(dofpe_p,1); %Fuerzas internas asociadas al flujo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fflux2 = zeros(dofpe_p,1); %Fuerzas internas asociadas al flujo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% if isfield(e_VG,'conshypMacro') && conshyp == 16
% %     fMul_Lag = zeros(dofpe_p,1); %Fuerzas internas asociadas al multiplicador de Lagrange
%     Lalfa = zeros(dofpe_p,dofpe_p);
% end

sigmaE_new = zeros(ntens,nPG); %AA
sigmaT_new = zeros(ntens,nPG); %AA
eps_new = zeros(ntens,nPG); % Corresponde a PG del elemento en la MICRO
eps_fluct = zeros(ntens,nPG); % Corresponde a PG del elemento en la MICRO
velflu_new = zeros(2,nPG);  % Corresponde a PG del elemento en la MICRO
velfluMicro_new = zeros(2,nPG);  % Corresponde a PG del elemento en la MICRO
mfluY_new = zeros(2,nPG);  % Corresponde a PG del elemento en la MICRO
mflu_new = zeros(1,nPG); % Corresponde a PG del elemento en la MICRO
hvar_new      = f_DefinicionhVar(conshyp,sihvarpg,nPG,e_VG.protype); %AA: add case 14
%AA22: Agrego para pasar de la MACRO a la MICRO las poropresiones y los
%gradientes de estas
phi_new = zeros(2,nPG); % Corresponde a PG del elemento en la MICRO
phi_fluct = zeros(2,nPG); % Corresponde a PG del elemento en la MICRO
p_fluct = zeros(4,1); % Corresponde a los nodos del elemento en la MICRO
p_M = zeros(1,nPG);
%AA22: Agrego para pasar de la MACRO a la MICRO las poropresiones y los
%gradientes de estas

% if esImplex
%     m_TensorTang = zeros(ntens,ntens,2*nPG);
% else
    m_TensorTang = zeros(ntens,ntens,nPG);
    m_Csig_eps= zeros(ntens,ntens,nPG);
    m_bw_eps= zeros(1,ntens,nPG);
    m_Cw_eps3= zeros(2,ntens,nPG);
    m_bsig_p= zeros(ntens,1,nPG);
    m_Mww_p= zeros(1,1,nPG);
    m_bw_p3=zeros(2,1,nPG);
    m_ksig_phi= zeros(ntens,2,nPG);
    m_kww_phi= zeros(1,2,nPG);
    m_kw_phi3= zeros(2,2,nPG);
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m_Cw_eps= zeros(2,ntens,nPG);
    m_bw_p=zeros(2,1,nPG);
    m_kw_phi= zeros(2,2,nPG);
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% end

hvar_old = reshape(hvar_old,sihvarpg,[]);
aux_var = reshape(aux_var,siavarpg,nPG);

%En el peso de gauss ya viene multiplicado el espesor.
m_pesoPG_d = m_DetJe_d.*wg;
m_pesoPG_p = m_DetJe_p.*wg;

ud = u(pos_d); % En la MICRO desplazamientos fluctuantes y en la MACRO desplazamientos
up = u(pos_p); % En la MICRO poropresiones fluctuantes y en la MACRO poropresiones

if isfield(e_VG,'conshypMacro') && conshyp == 15
    alfa = e_DatMatSet.alfa_penalidad;  
elseif isfield(e_VG,'conshypMacro') && conshyp == 16
    fMul_Lag = zeros(dofpe_p,1); %Fuerzas internas asociadas al multiplicador de Lagrange
    pos_lambda = e_DatElemSet.pos_lambda;
    ulambda = u(pos_lambda);
end

for iPG = 1:nPG

    e_VG.iPG = iPG;
    B = m_Be_d(:,:,iPG);
    DerivN = m_Dercae_p(:,:,iPG);
    N4 = m_FF_p(:,:,iPG);

    eps_fluct(:,iPG) = B*ud;
    phi_fluct(:,iPG) = DerivN*up;
    p_fluct(:,1) = up;

    % Deformacion aplicada a cada PG
    eps_new(:,iPG) = DefMacro + eps_fluct(:,iPG);
    % Gradiente de poropresiones aplicada a cada PG
    phi_new(:,iPG) = GradPorMacro + phi_fluct(:,iPG); 
    % Poropresiones aplicada a cada PG
    pNodo_new = PorMacro + (GradPorMacro'*(N4*(yy(1:2,1:4))')')  + N4*p_fluct;


    % Modelo constitutivo
    switch conshyp
        case {14,15,16} %Modelo bifasico con solido elastico (saturado) con metodo de la penalidad
            % Poropresiones aplicada a cada PG
%             pNodo_new = PorMacro + (GradPorMacro'*(N4*(yy(1:2,1:4))')')  + N4*p_fluct;
            [ct,BiotM,beta,PermK,sigmaE_new(:,iPG),sigmaT_new(:,iPG),mflu_new(1,iPG),velfluMicro_new(:,iPG)] = ...
            f_Rmap_Bif_Elast_Monoscale(eps_new(:,iPG),phi_new(:,iPG),pNodo_new,N4,e_DatMatSet);  %AA: cree funcion
            p_M(1,iPG) = pNodo_new; %En la MICRO no tiene sentido solo lo hago para UNIFORMIZAR LA SALIDA ?!?
        case {50,55}  %Modelo multiescala medio poroso  
            % Poropresiones en cada nodo del elemento
%             pNodo_new = PorMacro + (GradPorMacro'*(N4*(yy(1:2,1:4))')')  + N4*p_fluct;
            p_M(1,iPG) = pNodo_new;
            [e_TanOp,sigmaE_new(:,iPG),sigmaT_new(:,iPG),velflu_new(:,iPG),mflu_new(:,iPG),mfluY_new(:,iPG),hvar_new(:,iPG)] = ...
            f_RMap_ME_Bifase(eps_new(:,iPG),phi_new(:,iPG),p_M(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
        otherwise
            error('Matrices Elementales bifasicas: Modelo constitutivo no definido.')
    end
    if conshyp==14 
        %Se almacena para los tensor tangente constitutivos para realizar homogeneizacion y analisis de
        %bifurcacion
        if esImplex
            % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
            %tangentes implicitas para el analisis de bifurcacion.
            %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando esta seleccionado
            %el implex.
            m_TensorTang(:,:,iPG) = ct.Implex;
            %Se almacena los tensores tangentes constitutivo implicitos para analisis de bifurcacion como si fuera
            %PG adicionales, tantos como nPG. Se almacena en los ondices (:,:,nPG+1:2*nPG).
            %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
            %tercera dimension de m_TensorTang (size(m_TensorTang,3)).
            m_TensorTang(:,:,iPG+nPG) = ct.Impli;
            %En los calculos para el ensamblaje se utiliza el implex.
            ct = ct.Implex;
        else
            if isfield(ct,'Impli')
                ct = ct.Impli;
            end
        m_TensorTang(:,:,iPG) = ct;           
        end
        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*ct*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + B'*BiotM*N4*m_pesoPG_p(iPG); %Matriz de acople
        Qps = Qps + N4'*BiotM'*B*m_pesoPG_p(iPG); %Matriz de acople %Qps = Qps + N4*BiotM'*B*m_pesoPG_p(iPG);
        Sww = Sww + beta*(N4'*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); % Matriz de permeabilidad %Hww = Hww + DerivN*PermK*DerivN'*m_pesoPG_p(iPG);
        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); % Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velfluMicro_new(:,iPG)*m_pesoPG_p(iPG); % Fuerzas internas debidas a la velocidad de filtracion del fluido

        velflu_new(:,iPG) = velfluMicro_new(:,iPG);
        
    elseif conshyp==15 || conshyp==16
        %Se almacena para los tensor tangente constitutivos para realizar homogeneizacion y analisis de
        %bifurcacion
        if esImplex
            % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
            %tangentes implicitas para el analisis de bifurcacion.
            %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando esta seleccionado
            %el implex.
            m_TensorTang(:,:,iPG) = ct.Implex;
            %Se almacena los tensores tangentes constitutivo implicitos para analisis de bifurcacion como si fuera
            %PG adicionales, tantos como nPG. Se almacena en los ondices (:,:,nPG+1:2*nPG).
            %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
            %tercera dimension de m_TensorTang (size(m_TensorTang,3)).
            m_TensorTang(:,:,iPG+nPG) = ct.Impli;
            %En los calculos para el ensamblaje se utiliza el implex.
            ct = ct.Implex;
        else
            if isfield(ct,'Impli')
                ct = ct.Impli;
            end
        m_TensorTang(:,:,iPG) = ct;           
        end
        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*ct*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + B'*BiotM*N4*m_pesoPG_p(iPG); %Matriz de acople
        Qps = Qps + N4'*BiotM'*B*m_pesoPG_p(iPG); %Matriz de acople %Qps = Qps + N4*BiotM'*B*m_pesoPG_p(iPG);
        Sww = Sww + beta*(N4'*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); % Matriz de permeabilidad %Hww = Hww + DerivN*PermK*DerivN'*m_pesoPG_p(iPG);
        Lalfa = Lalfa + N4'*N4*m_pesoPG_p(iPG);
        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); % Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velfluMicro_new(:,iPG)*m_pesoPG_p(iPG); % Fuerzas internas debidas a la velocidad de filtracion del fluido

        if isfield(e_VG,'conshypMacro')
            velflu_new(:,iPG) = velfluMicro_new(:,iPG);%-(mflu_new(:,iPG)*N4*(yy(1:2,1:4))')';
            mfluY_new(:,iPG) = -(mflu_new(:,iPG)*N4*(yy(1:2,1:4))')';
        else
            velflu_new(:,iPG) = velfluMicro_new(:,iPG);
        end
%     elseif conshyp==16
%         %Se almacena para los tensor tangente constitutivos para realizar homogeneizacion y analisis de
%         %bifurcacion
%         if esImplex
%             % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
%             %tangentes implicitas para el analisis de bifurcacion.
%             %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando esta seleccionado
%             %el implex.
%             m_TensorTang(:,:,iPG) = ct.Implex;
%             %Se almacena los tensores tangentes constitutivo implicitos para analisis de bifurcacion como si fuera
%             %PG adicionales, tantos como nPG. Se almacena en los ondices (:,:,nPG+1:2*nPG).
%             %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
%             %tercera dimension de m_TensorTang (size(m_TensorTang,3)).
%             m_TensorTang(:,:,iPG+nPG) = ct.Impli;
%             %En los calculos para el ensamblaje se utiliza el implex.
%             ct = ct.Implex;
%         else
%             if isfield(ct,'Impli')
%                 ct = ct.Impli;
%             end
%         m_TensorTang(:,:,iPG) = ct;           
%         end
%         % Calculo de las sumbatrices del la matriz de rigidez acoplada
%         Kss = Kss + B'*ct*B*m_pesoPG_d(iPG); %Matriz de rigidez
%         Qsp = Qsp + B'*BiotM*N4*m_pesoPG_p(iPG); %Matriz de acople
%         Qps = Qps + N4'*BiotM'*B*m_pesoPG_p(iPG); %Matriz de acople %Qps = Qps + N4*BiotM'*B*m_pesoPG_p(iPG);
%         Sww = Sww + beta*(N4'*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad
%         Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); % Matriz de permeabilidad %Hww = Hww + DerivN*PermK*DerivN'*m_pesoPG_p(iPG);
%         % Calculo de fuerzas internas (P/residuo)
%         fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); % Fuerzas internas debidas a las tensiones totales
%         fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas al contenido de masa del fluido
%         fflux_V = fflux_V + DerivN'*velfluMicro_new(:,iPG)*m_pesoPG_p(iPG); % Fuerzas internas debidas a la velocidad de filtracion del fluido
% 
%         velflu_new(:,iPG) = velfluMicro_new(:,iPG);
    elseif conshyp==50
%         m_TensorTang(:,:,iPG)= e_TanOp.m_Csig_eps;
        m_Csig_eps(:,:,iPG)= e_TanOp.m_Csig_eps;
        m_bw_eps(:,:,iPG)= e_TanOp.m_bw_eps;
        m_Cw_eps3(:,:,iPG)= e_TanOp.m_Cw_eps3;
        m_bsig_p(:,:,iPG)= e_TanOp.m_bsig_p;
        m_Mww_p(:,:,iPG)= e_TanOp.m_Mww_p;
        m_bw_p3(:,:,iPG)= e_TanOp.m_bw_p3;
        m_ksig_phi(:,:,iPG)= e_TanOp.m_ksig_phi;
        m_kww_phi(:,:,iPG)= e_TanOp.m_kww_phi;
        m_kw_phi3(:,:,iPG)= e_TanOp.m_kw_phi3;
        
        m_Cw_eps(:,:,iPG)= e_TanOp.m_Cw_eps;
        m_bw_p(:,:,iPG)= e_TanOp.m_bw_p;
        m_kw_phi(:,:,iPG)= e_TanOp.m_kw_phi;

        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*m_Csig_eps(:,:,iPG)*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + (B'*m_ksig_phi(:,:,iPG)*DerivN+B'*m_bsig_p(:,:,iPG)*N4)*m_pesoPG_p(iPG); %Matriz de acoplez
        Qps = Qps + (-DerivN'*m_Cw_eps3(:,:,iPG)*B+N4'*m_bw_eps(:,:,iPG)*B)*m_pesoPG_p(iPG); %Matriz de acople
        Kpp = Kpp + (-DerivN'*m_kw_phi3(:,:,iPG)*DerivN-DerivN'*m_bw_p3(:,:,iPG)*N4+...
        N4'*m_kww_phi(:,:,iPG)*DerivN+N4'*m_Mww_p(:,:,iPG)*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad + Matriz de permeabilidad!!!!?????       

        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); %Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); %Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velflu_new(:,iPG)*m_pesoPG_p(iPG); %Fuerzas internas debidas a la velocidad de filtracion del fluido
        fflux_mY = fflux_mY + DerivN'*mfluY_new(:,iPG)*m_pesoPG_p(iPG); %Fuerzas internas debidas a la velocidad de filtracion del fluido
    end %if(conshyp)   
end %for(ipg)

if conshyp==14
%     Kpp = Sww + (theta*Dtime*(Hww+alfa*Lalfa));
    Kpp = Sww + theta*Dtime*Hww;
%     fflux  = fflux - fflux_m - ((theta*Dtime)*(fflux_V+(alfa*Lalfa*up))); %Fuerzas internas asociadas al flujo
    fflux  = fflux - fflux_m - (theta*Dtime)*fflux_V; %Fuerzas internas asociadas al flujo
%     fflux  = fflux - fflux_m - (theta*Dtime)*fflux_V; %Fuerzas internas asociadas al flujo
%     fflux  = fflux - Qps*ud - (Sww+theta*Dtime*Hww)*up; %Fuerzas internas asociadas al flujo
%     fstreT = fstreT + Kss*ud - Qsp*up;
    % Vector de fuerzas elemental acoplado
    Fe = [fstreT ; fflux ; zeros(4,1)];

    % Matriz de rigidez elemental acoplada
    Ke = [Kss -Qsp zeros(16,4);
    -Qps -Kpp zeros(4,4) ;
    zeros(4,20) eye(4,4)];
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand(Ke,Fe,dofpe) ;

    % Se ordena las matrices como vectores columnas
    eps_new    = eps_new(:);
    eps_fluct  = eps_fluct(:);
%     e_TanOp = struct('m_TensorTang',m_TensorTang);
elseif conshyp==15
%     Kpp = Sww + (theta*Dtime*(Hww+alfa*Lalfa));
    Kpp = Sww + Dtime*alfa*Lalfa + theta*Dtime*Hww;
%     fflux  = fflux - fflux_m - ((theta*Dtime)*(fflux_V+(alfa*Lalfa*up))); %Fuerzas internas asociadas al flujo
    fflux  = fflux - fflux_m - Dtime*alfa*Lalfa*up - (theta*Dtime)*fflux_V; %Fuerzas internas asociadas al flujo
%     fflux  = fflux - fflux_m - (theta*Dtime)*fflux_V; %Fuerzas internas asociadas al flujo
%     fflux  = fflux - Qps*ud - (Sww+theta*Dtime*Hww)*up; %Fuerzas internas asociadas al flujo
%     fstreT = fstreT + Kss*ud - Qsp*up;
    % Vector de fuerzas elemental acoplado
    Fe = [fstreT ; fflux ; zeros(4,1)];

    % Matriz de rigidez elemental acoplada
    Ke = [Kss -Qsp zeros(16,4);
    -Qps -Kpp zeros(4,4) ;
    zeros(4,20) eye(4,4)];
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand(Ke,Fe,dofpe) ;

    % Se ordena las matrices como vectores columnas
    eps_new    = eps_new(:);
    eps_fluct  = eps_fluct(:);
%     e_TanOp = struct('m_TensorTang',m_TensorTang);
elseif conshyp==16
%     Kpp = Sww + (theta*Dtime*(Hww+alfa*Lalfa));
    Kpp = Sww + theta*Dtime*Hww;
%     Lalfa = theta*Dtime*Lalfa;      
    fflux  = fflux - fflux_m - (theta*Dtime)*fflux_V - Lalfa*ulambda; %Fuerzas internas asociadas al flujo
    fMul_Lag = fMul_Lag - Lalfa*up; %Fuerzas internas asociadas al multiplicador de Lagrange
    % Vector de fuerzas elemental acoplado
    Fe = [fstreT ; fflux ; zeros(4,1); fMul_Lag; zeros(4,1)];

    % Matriz de rigidez elemental acoplada
%     Ke = [Kss -Qsp zeros(16,12);
%     -Qps -Kpp zeros(4,4) -Lalfa zeros(4,4);
%     zeros(4,20) eye(4,4) zeros(4,8);
%     zeros(4,16) -Lalfa zeros(4,4) 0.0000000001*eye(4,4) zeros(4,4);
%     zeros(4,28) 0.0000000001*eye(4,4)];
    Ke = [Kss -Qsp zeros(16,12);
    -Qps -Kpp zeros(4,4) -Lalfa zeros(4,4);
    zeros(4,20) eye(4,4) zeros(4,8);
    zeros(4,16) -Lalfa zeros(4,4) 0.00000000000001*eye(4,4) zeros(4,4);
    zeros(4,28) 0.00000000000001*eye(4,4)];
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand_ML(Ke,Fe,dofpe) ;

    % Se ordena las matrices como vectores columnas
    eps_new    = eps_new(:);
    eps_fluct  = eps_fluct(:);
%     e_TanOp = struct('m_TensorTang',m_TensorTang);
elseif conshyp==50
    velflu_new = velflu_new - mfluY_new; % VELOCIDAD DE FILTRACION HOMOGENEIZADO PARA SALIDA DE DATOS.
    fflux = fflux - Qps*ud - Kpp*up; %Fuerzas internas asociadas al flujo
%         fflux  = fflux - fflux_m + fflux_V; %Fuerzas internas asociadas al flujo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fflux2  = fflux2 - fflux_m - ((theta*Dtime)*fflux_V + fflux_mY); %Fuerzas internas asociadas al flujo
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         fflux  = fflux - fflux_m - (theta*DTime)*fflux_V; %Fuerzas internas asociadas al flujo
    % Vector de fuerzas elemental acoplado
    Fe = [fstreT ; fflux ; zeros(4,1)];
    % Matriz de rigidez elemental acoplada
    % NO ENTIENDO PORQUE TENGO QUE PONER EN NEGATIVO LA SEGUNDA LINEA
    % A MI PARECER TODO DEBE SER POSITIVO PERO CON LOS SIGNOS COMO
    % ESTA CORRE, FALTA REVISAR SI LOS RESULTADOS SON CORRECTOS!!!!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % PODRIA SER QUE EL NEGATIVO ES NECESARIO PORQUE EN LA MICRO LO TOMO NEGATIVO,
    % ENTONCES AL RESOLVER LAS VARIABLES FLUCTUANTES ESTOY "INTRODUCIENDO" ESE NEGATIVO,
    % QUE SE TRASLADA A LOS TENSORES HOMOGENEIZADOS MACRO QUE "CONFORMAN" LA 
    % LINEA SEGUNDA LINEA. ENTONCES, AHORA EN LA MACRO ESTOS DEBERIAN SER POSITIVOS
    % LO CUAL CORRIJO CON EL SIGNO MENOS?
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Ke = [Kss  Qsp zeros(16,4);
         -Qps -Kpp zeros(4,4) ;
         zeros(4,20) eye(4,4)];
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand(Ke,Fe,dofpe) ;
end %if(conshyp) 

% Se ordena las matrices como vectores columnas
eps_new    = eps_new(:); 
eps_fluct  = eps_fluct(:); 
sigmaE_new  = sigmaE_new(:); %AA (por PG)
sigmaT_new  = sigmaT_new(:); %AA (por PG)
mflu_new  = mflu_new(:); %AA (por PG)
velflu_new  = velflu_new(:); %AA (por PG)
mfluY_new = mfluY_new(:);
phi_new    = phi_new(:); 
phi_fluct  = phi_fluct(:); 
pNodo_new    = pNodo_new(:); 
p_fluct  = p_fluct(:); 
p_M = p_M(:);
hvar_new   = hvar_new(:);
aux_var    = aux_var(:);

