function [kt,fint,eps_new,phi_new,p_new,eps_fluct,phi_fluct,p_fluct,...
            sigmaE_new,sigmaT_new,mflu_new,velflu_sta,velflu_total,...
            velflu_new,hvar_new,aux_var,m_sig_eps,m_chi_eps,m_V_eps,...
            m_sig_p,m_chi_p,m_V_p,m_sig_phi,m_chi_phi,m_V_phi]=...
            f_MatElem_BifaseMulTSc(u,u_old,yy,hvar_old,aux_var,...
            m_Be_d,m_DetJe_d,m_Dercae_p,m_DetJe_p,m_FF_d,m_FF_p,DefMacro,...
            GradPorMacro_dup,GradPorMacro_up_n,PorMacro,...
            DefMacro_new,GradPorMacro_up_new,PorMacro_new,...
            e_DatElemSet,e_DatMatSet,e_VG)
        
%AA: Cree funcion
% Variable globales
ntens = e_VG.ntens; %Numero de tensiones
theta = e_VG.theta; %Factor de regla de punto medio integracion de Euler
Dtime = e_VG.Dtime; %Paso de tiempo 
dofpe = e_DatElemSet.dofpe; %Grados de libertad del elemento 
dofpe_d = e_DatElemSet.dofpe_d; %Grados de libertad en desplazamientos
dofpe_p = e_DatElemSet.dofpe_p; %Grados de libertad en poro presiones
nPG = e_DatElemSet.npg; %Numero de puntos de Gauss
wg = e_DatElemSet.wg;
pos_d = e_DatElemSet.pos_d;
pos_p = e_DatElemSet.pos_p;

sihvarpg = e_DatMatSet.sihvarpg;
siavarpg = e_DatMatSet.siavarpg;
conshyp = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;

% Inicializaciones
%Matriz de rigidez
Kss = zeros(dofpe_d,dofpe_d);
%Matriz de acople
Qsp = zeros(dofpe_d,dofpe_p);
%Matriz de acople traspuesta
Qps = zeros(dofpe_p,dofpe_d);
%Matriz de compresibilidad
Sww = zeros(dofpe_p,dofpe_p);
%Matriz de permeabilidad
Hww = zeros(dofpe_p,dofpe_p);
%Matriz de 'rigidez' para poro presiones. Combina Sww y Hww
Kpp = zeros(dofpe_p,dofpe_p);

%Fuerzas internas debidas a las tensiones totales
fstreT = zeros(dofpe_d,1); 
%'Fuerzas' internas asociadas al flujo por el contenido de masa del fluido
fflux_m = zeros(dofpe_p,1);
%'Fuerzas' internas asociadas al flujo por el termino estatico de la velocidad de filtracion
fflux_V = zeros(dofpe_p,1); 
% % % %Fuerzas internas asociadas al flujo por el termino dinamico de la velocidad de filtracion
% % % fflux_mY = zeros(dofpe_p,1);
%Fuerzas internas asociadas al flujo. Combinacion de las anteriores
fflux = zeros(dofpe_p,1); 

%############################################################################################
%VARIABLES PRIMITIVAS
%Delta de deformaciones
eps_new = zeros(ntens,nPG); 
%Deformaciones  al paso "n+1"
eps_n1 = zeros(ntens,nPG); 
%Delta del gradiente de poro presiones
phi_new = zeros(2,nPG);
%Gradiente de poro presiones al paso "n"
phi_new_up_n = zeros(2,nPG);
%Gradiente de poro presiones al paso "n+1"
phi_new_up_n1 = zeros(2,nPG);
%Delta de poro presiones
p_new = zeros(1,nPG);
%Poro presiones al paso "n+1"
p_n1 = zeros(1,nPG);


%VARIABLES PRIMITIVAS FLUCTUANTES
%Delta de micro-deformaciones fluctuantes o delta de macro-deformaciones al paso "n+1"
eps_fluct = zeros(ntens,nPG);
%Micro-deformaciones fluctuantes o delta de macro-deformaciones
eps_fluct_n1 = zeros(ntens,nPG);
%Delta de gradiente de micro-poro presiones fluctuantes o delta de gradiente de macro-poro presiones
phi_fluct = zeros(2,nPG);
%Gradiente de micro-poro presiones fluctuantes al paso "n" o gradiente de
%macro-poro presiones al paso "n"
phi_fluct_up_n = zeros(2,nPG);
%Gradiente de micro-poro presiones fluctuantes al paso "n+1" o gradiente de
%macro-poro presiones al paso "n+1"
phi_fluct_up_n1 = zeros(2,nPG);
%Delta de micro-poro presiones fluctuantes o delta de macro-poro presiones
p_fluct = zeros(4,1);
%Micro-poro presiones fluctuantes o delta de macro-poro presiones al paso "n+1"
p_fluct_n1 = zeros(4,1);
%############################################################################################

%############################################################################################
%VARIABLES DUALES
%Delta de tensiones efectivas
sigmaE_new = zeros(ntens,nPG);
%Delta de tensiones totales
sigmaT_new = zeros(ntens,nPG);
%Delta de tasa del contenido de masa del fluido
mflu_new = zeros(1,nPG); 
%Componente estatica de la velocidad de filtracion  al tiempo "n+theta"
velflu_sta = zeros(2,nPG);
%Delta de la componente dinamica de la velocidad de filtracion
velflu_dyn = zeros(2,nPG); 
%Velocidad de filtracion  al tiempo "n+theta"
velflu_new = zeros(2,nPG); 

%Velocidad de filtracion estatica al tiempo "n+1"
velflu_n1 = zeros(2,nPG); 
%Velocidad de filtracion dinamica al tiempo "n+1"
mflu_n1 = zeros(1,nPG);

%Velocidad de filtracion  al tiempo "n+1"
velflu_total = zeros(2,nPG); 


%############################################################################################
%Agregue Fext
%Variables historicas de cada modelo constitutivo
hvar_new = f_DefinicionhVar(conshyp,sihvarpg,nPG,e_VG.protype); %AA: add case 14
%############################################################################################
% if esImplex
%     m_TensorTang = zeros(ntens,ntens,2*nPG);
% else
%############################################################################################
%OPERADORES TANGENTES
%Operador tangente constitutivo dsigma/deps
    m_sig_eps= zeros(ntens,ntens,nPG);
%Operador tangente constitutivo dchi/deps
    m_chi_eps= zeros(1,ntens,nPG);
%Operador tangente constitutivo dV/deps
    m_V_eps= zeros(2,ntens,nPG);
%Operador tangenteconstitutivo dsigma/dp
    m_sig_p= zeros(ntens,1,nPG);
%Operador tangente constitutivo dchi/dp    
    m_chi_p= zeros(1,1,nPG);
%Operador tangente constitutivo dV/dp       
    m_V_p=zeros(2,1,nPG);
%Operador tangente constitutivo dsigma/dphi    
    m_sig_phi= zeros(ntens,2,nPG);
%Operador tangente constitutivo dchi/dphi    
    m_chi_phi= zeros(1,2,nPG);
%Operador tangente constitutivo dV/dphi    
    m_V_phi= zeros(2,2,nPG);
%end
%############################################################################################

%############################################################################################
%VARIABLES DEPENDIENTES DEL MODELO CONSTITUTIVO
%Estructura de multiples variables dependiente del modelo constitutivo
hvar_old = reshape(hvar_old,sihvarpg,[]);
aux_var = reshape(aux_var,siavarpg,nPG);
%############################################################################################

%############################################################################################
%PESO DE INTEGRACION
%En el peso de gauss ya viene multiplicado el espesor.
m_pesoPG_d = m_DetJe_d.*wg;
m_pesoPG_p = m_DetJe_p.*wg;
%############################################################################################

%############################################################################################
%VARIABLES PRIMALES
%Delta de VARIABLES PRIMALES
du = u - u_old;
%En la MICRO incognitas fluctuantes y en la MACRO incognitas
%Delta de desplazamientos (u_n+1 - u_n)
dud = du(pos_d);
%Delta de poro presiones (p_n+1 - p_n)
dup = du(pos_p);
%Poro presiones al paso de tiempo "n"
ud_old = u_old(pos_d); % d_old    
%Poro presiones al paso de tiempo "n+1"
ud_new = u(pos_d); % d_new   
%Poro presiones al paso de tiempo "n"
up_old = u_old(pos_p); % p_old    
%Poro presiones al paso de tiempo "n+1"
up_new = u(pos_p); % p_new   
%############################################################################################

%############################################################################################
%CONSTANTES Y MATRICES NECESARIAS S/MODELO CONSTITUTIVO
if conshyp == 14 %Caso conshyp=14 MONO-ESCALA o RVE. 
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %OJO!! MONO-ESCALA debe ser beta_factor = 1 siempre 
    %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %Factor binario segun tipo de expansion del argumento en la ecuacion
    %constituva Full Order Expanded = 1 y Zero Order Expanded = 0
    beta_factor1 = 1; %Por ahora solo puedo modificarlo desde aca 
                     %y no desde el archivo de entrada
    beta_factor2 = 1; %Por ahora solo puedo modificarlo desde aca 
    %y no desde el archivo de entrada
% elseif conshyp == 15 % caso conshyp=15 RVE de caso MULTIESCALA
%     %Inicializacion
%     %Matriz que contiene factor de penalidad
%     Lalfa = zeros(dofpe_p,dofpe_p); 
%     %Coeficiente de penalidad
%     alfa = e_DatMatSet.alfa_penalidad;
%     %Factor binario segun tipo de expansion del argumento en la ecuacion
%     %constituva Full Order Expanded = 1 y Zero Order Expanded = 0
%     beta_factor = e_DatMatSet.beta_factor;
elseif conshyp == 15 % caso conshyp=15 RVE de caso MULTIESCALA    
    %Inicializacion
    Llambda_p = zeros(dofpe_p,1); 
    %Fuerzas internas asociadas al multiplicador de Lagrange por
    %poro presiones
    fMul_Lag_p = 0;
    %Posicion de los terminos asociados a multiplicador de Lagrange
    pos_lambda_p = e_DatElemSet.pos_lambda_p;
    %Vector de terminos de delta de multiplicadores de Lagrange por
    %micro-poro presiones fluctuantes
    dulambda_p = du(pos_lambda_p);
    %Factor binario segun tipo de expansion del argumento en la ecuacion
    %constituva Full Order Expanded = 1 y Zero Order Expanded = 0
    beta_factor1 = e_DatMatSet.beta_factor1;
    beta_factor2 = e_DatMatSet.beta_factor2;
elseif conshyp == 16% caso conshyp=16 RVE de caso MULTIESCALA
    %Inicializacion
   %Vector que contiene terminos adicionales por el Mult. Lagrange por
    %desplazamientos
    Llambda_u = zeros(dofpe_d,1); 
    %Fuerzas internas asociadas al multiplicador de Lagrange por
    %desplazamientos
    fMul_Lag_u = zeros(2,1);
    %Posicion de los terminos asociados a multiplicador de Lagrange por
    %desplazamientos
    pos_lambda_u = e_DatElemSet.pos_lambda_u;
    %Vector de terminos de delta de multiplicadores de Lagrange por
    %micro-desplazamientos fluctuantes
    dulambda_u = du(pos_lambda_u);
    %Factor binario segun tipo de expansion del argumento en la ecuacion
    %constituva Full Order Expanded = 1 y Zero Order Expanded = 0
    beta_factor1 = e_DatMatSet.beta_factor1;
    beta_factor2 = e_DatMatSet.beta_factor2;
elseif conshyp == 17 % caso conshyp=17 RVE de caso MULTIESCALA    
    %Inicializacion
    %Vector que contiene terminos adicionales por el Mult. Lagrange por
    %desplazamientos
    Llambda_u = zeros(dofpe_d,1); 
    %Vector que contiene terminos adicionales por el Mult. Lagrange por
    %poro presiones
    Llambda_p = zeros(dofpe_p,1); 
    %Fuerzas internas asociadas al multiplicador de Lagrange por
    %desplazamientos
    fMul_Lag_u = zeros(2,1);
    %Fuerzas internas asociadas al multiplicador de Lagrange por
    %poro presiones
    fMul_Lag_p = 0;
    %Posicion de los terminos asociados a multiplicador de Lagrange por
    %desplazamientos
    pos_lambda_u = e_DatElemSet.pos_lambda_u;
    %Posicion de los terminos asociados a multiplicador de Lagrange
    pos_lambda_p = e_DatElemSet.pos_lambda_p;
    %Vector de terminos de multiplicadores de Lagrange por
    %micro-desplazamientos fluctuantes al paso de tiempo "n"
%     ulambda_u = u_old(pos_lambda_u);
    %Vector de terminos de delta de multiplicadores de Lagrange por
    %micro-desplazamientos fluctuantes
    dulambda_u = du(pos_lambda_u);
    %Vector de terminos de multiplicadores de Lagrange por
    %micro-poro presiones fluctuantes al paso de tiempo "n"
%     ulambda_p = u_old(pos_lambda_p);
    %Vector de terminos de delta de multiplicadores de Lagrange por
    %micro-poro presiones fluctuantes
    dulambda_p = du(pos_lambda_p);
    %Factor binario segun tipo de expansion del argumento en la ecuacion
    %constituva Full Order Expanded = 1 y Zero Order Expanded = 0
    beta_factor1 = e_DatMatSet.beta_factor1;
    beta_factor2 = e_DatMatSet.beta_factor2;
else % caso conshyp=50 MACRO-ESCALA de caso MULTIESCALA
    %Macro-escala. Siempre debe ser igual 1
    beta_factor1 = 1;
    beta_factor2 = 1;
end
%############################################################################################
%############################################################################################
%BUCLE DE INTEGRACION EN PUNTOS DE GAUSS
for iPG = 1:nPG
    %Punto de integracion en estudio
    e_VG.iPG = iPG;
    %Matriz desplazamiento-deformacion
    B = m_Be_d(:,:,iPG);
    %Matriz que relacion gradiente de poro presiones con poro presiones
    DerivN = m_Dercae_p(:,:,iPG);
    %Matriz de funcion de interpolacion de elemento cuadrilatero
    %bicuadratico
    N8 = m_FF_d(:,:,iPG);
    %Matriz de funcion de interpolacion de elemento cuadrilatero bilineal
    N4 = m_FF_p(:,:,iPG);
        
    %VARIABLES PRIMITIVAS FLUCTUANTES    
    %Delta de micro-deformaciones fluctuantes o delta de macro-deformaciones
    eps_fluct(:,iPG) = B*dud;
    %Micro-deformaciones fluctuantes o delta de macro-deformaciones
    eps_fluct_n1(:,iPG) = B*ud_new;
    %Delta de gradiente de micro-poro presiones fluctuantes o delta de gradiente de macro-poro presiones
    phi_fluct(:,iPG) = DerivN*dup;
    %Gradiente de micro-poro presiones fluctuantes al paso "n" o
    %gradiente de macro-poro presiones al paso "n"
    phi_fluct_up_n(:,iPG) = DerivN*up_old;
    %Gradiente de micro-poro presiones fluctuantes al paso "n+1" o gradiente de
    %macro-poro presiones al paso "n+1"
    phi_fluct_up_n1(:,iPG) = DerivN*up_new;
    %Delta de micro-poro presiones fluctuantes o delta de macro-poro presiones
    p_fluct(:,1) = dup;
    %Micro-poro presiones fluctuantes o delta de macro-poro presiones al paso "n+1"
    p_fluct_n1(:,1) = up_new;

    %VARIABLES PRIMITIVAS
    %Delta de deformaciones
    eps_new(:,iPG) = DefMacro + eps_fluct(:,iPG);
    %Deformaciones  al paso "n+1"
    eps_n1(:,iPG) =  DefMacro_new + eps_fluct_n1(:,iPG);
    %Delta del gradiente de poro presiones "
    phi_new(:,iPG) = GradPorMacro_dup + phi_fluct(:,iPG);      
    %Gradiente de poro presiones al paso "n"
    phi_new_up_n(:,iPG) = GradPorMacro_up_n + phi_fluct_up_n(:,iPG); 
    %Gradiente de micro-poro presiones fluctuantes al paso "n+1" o gradiente de
    %macro-poro presiones al paso "n+1"
    phi_new_up_n1(:,iPG) =  GradPorMacro_up_new + phi_fluct_up_n1(:,iPG); 
    %Micro: depende del modelo constitutivo utilizado. Si el argumento es
    %Full Order Expanded (beta=1) o Zero Order Expanded (beta=0)
    p_new(1,iPG) = PorMacro + beta_factor1*((GradPorMacro_dup'*(N4*(yy(1:2,1:4))')'))  + beta_factor2*(N4*p_fluct);
    %Micro-poro presiones fluctuantes o delta de macro-poro presiones al
    %paso "n+1".Si el argumento es
    %Full Order Expanded (beta=1) o Zero Order Expanded (beta=0)
    p_n1(:,iPG) = PorMacro_new + beta_factor1*((GradPorMacro_up_new'*(N4*(yy(1:2,1:4))')'))  + beta_factor2*(N4*p_fluct_n1);

%############################################################################################
    % Modelos constitutivos
    switch conshyp
        case {14} %Modelo bifasico con solido elastico (saturado) - Caso Mono-escala y para micro-escala en caso multi-escala
            %No permitiria diferenciar la p para chi de la p para sigma
            %como esta propuesto ahora
            p_new_sigma=p_new;
            
            [ct,BiotM,beta,PermK,sigmaE_new(:,iPG),sigmaT_new(:,iPG),...
                mflu_new(1,iPG),velflu_sta(:,iPG),mflu_n1(:,iPG),velflu_n1(:,iPG)] = ...
            f_Rmap_Bif_Elast_Monoscale(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),p_new_sigma(:,iPG),eps_n1(:,iPG),phi_new_up_n1(:,iPG),...
            p_n1(:,iPG),e_DatMatSet,theta); 
        
        %###############################################################################
        case {15,16,17} %Modelo bifasico con solido elastico (saturado) con metodo de la penalidad o multiplicadores de Lagrange
            %Caso Multi-escala
           p_new_sigma=p_new;
           
           [ct,BiotM,beta,PermK,sigmaE_new(:,iPG),sigmaT_new(:,iPG),...
               mflu_new(1,iPG),velflu_sta(:,iPG),mflu_n1(:,iPG),velflu_n1(:,iPG)] = ...
            f_Rmap_Bif_ElastME(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),p_new_sigma(:,iPG),eps_n1(:,iPG),phi_new_up_n1(:,iPG),...
            p_n1(:,iPG),e_DatMatSet,theta); 
        
        %###############################################################################
        case {50,55}  %Modelo multiescala medio poroso saturado  
            [e_TanOp,sigmaE_new(:,iPG),sigmaT_new(:,iPG),mflu_new(:,iPG),...
                velflu_sta(:,iPG),velflu_total(:,iPG),velflu_new(:,iPG),hvar_new(:,iPG)] = ...
            f_RMap_ME_Bifase(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),eps_n1(:,iPG),phi_new_up_n1(:,iPG),...
            p_n1(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
        
        case 60  %Modelo multiescala medio poroso saturado - Resolucion analitica 
            %modelo de Taylor completo
            [e_TanOp,~,~,~,~,~,~,hvar_new(:,iPG)] = ...
            f_RMap_ME_Bifase(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),eps_n1(:,iPG),phi_new_up_n1(:,iPG),...
            p_n1(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
        
            p_new_sigma=p_new;
           
            [ct,BiotM,beta,PermK,sigmaE_new(:,iPG),sigmaT_new(:,iPG),...
               mflu_new(1,iPG),velflu_new(:,iPG),velflu_total(:,iPG)] = ...
            f_Rmap_Analitico_TayComp(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),p_new_sigma(:,iPG),phi_new_up_n1(:,iPG),e_DatMatSet,theta,Dtime); 
        
        case 61  %Modelo multiescala medio poroso saturado - Resolucion analitica 
            %modelo Periodico en micro desplazamientos fluctuantes y Taylor 
            %en micro poro presiones fluctuantes
            [e_TanOp,~,~,~,~,~,~,hvar_new(:,iPG)] = ...
            f_RMap_ME_Bifase(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),eps_n1(:,iPG),phi_new_up_n1(:,iPG),...
            p_n1(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
        
            p_new_sigma=p_new;
           
            [ct,BiotM,beta,PermK,sigmaE_new(:,iPG),sigmaT_new(:,iPG),...
               mflu_new(1,iPG),velflu_new(:,iPG),velflu_total(:,iPG)] = ...
            f_Rmap_Analitico_PerTay(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),p_new_sigma(:,iPG),phi_new_up_n1(:,iPG),e_DatMatSet,theta,Dtime); 
        case 62  %Modelo multiescala medio poroso saturado - Resolucion analitica 
            %modelo  en Taylor micro-desplazamientos fluctuantes y Periodico 
            %en micro-poro presiones fluctuantes
            [e_TanOp,~,~,~,~,~,~,hvar_new(:,iPG)] = ...
            f_RMap_ME_Bifase(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),eps_n1(:,iPG),phi_new_up_n1(:,iPG),...
            p_n1(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
        
            p_new_sigma=p_new;
           
            [ct,BiotM,beta,PermK,sigmaE_new(:,iPG),sigmaT_new(:,iPG),...
               mflu_new(1,iPG),velflu_new(:,iPG),velflu_total(:,iPG)] = ...
            f_Rmap_Analitico_TayPer(eps_new(:,iPG),phi_new(:,iPG),phi_new_up_n(:,iPG),...
            p_new(:,iPG),p_new_sigma(:,iPG),phi_new_up_n1(:,iPG),e_DatMatSet,theta,Dtime); 
        
        otherwise
            error('Matrices Elementales bifasicas: Modelo constitutivo no definido.')
    end
    if conshyp==14 %Modelo bifasico con solido elastico (saturado) - Caso Mono-escala y para micro-escala en caso multi-escala
        m_sig_eps(:,:,iPG) = ct;  
        
        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*ct*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + B'*BiotM*N4*m_pesoPG_p(iPG); %Matriz de acople
        Qps = Qps + N4'*BiotM'*B*m_pesoPG_p(iPG); %Matriz de acople %Qps = Qps + N4*BiotM'*B*m_pesoPG_p(iPG);
        Sww = Sww + beta*(N4'*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); % Matriz de permeabilidad %Hww = Hww + DerivN*PermK*DerivN'*m_pesoPG_p(iPG);
       
        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); % Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velflu_sta(:,iPG)*m_pesoPG_p(iPG); % Fuerzas internas debidas a la velocidad de filtracion del fluido

        if isfield(e_VG,'conshypMacro')
             velflu_new(:,iPG) = velflu_sta(:,iPG)-((mflu_new(:,iPG)/Dtime)*N4*(yy(1:2,1:4))')';
             velflu_total(:,iPG) = velflu_n1(:,iPG)-((mflu_new(:,iPG)/Dtime)*N4*(yy(1:2,1:4))')';
        else
            velflu_new(:,iPG) = velflu_sta(:,iPG);
            velflu_total(:,iPG) = velflu_n1(:,iPG);
        end      
    elseif conshyp==15  %Modelo bifasico con solido elastico (saturado) con metodo de
        %multiplicadores de Lagrange en poropresiones. Caso Multi-escala
        m_sig_eps(:,:,iPG) = ct;   
        
        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*ct*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + B'*BiotM*N4*m_pesoPG_p(iPG); %Matriz de acople
        Qps = Qps + N4'*BiotM'*B*m_pesoPG_p(iPG); %Matriz de acople %Qps = Qps + N4*BiotM'*B*m_pesoPG_p(iPG);
        Sww = Sww + beta*(N4'*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); % Matriz de permeabilidad %Hww = Hww + DerivN*PermK*DerivN'*m_pesoPG_p(iPG);
        Llambda_p = Llambda_p + N4'*m_pesoPG_p(iPG); %Vector adicional para metodo de Mult. Lagrange por por opresiones
        
        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); % Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velflu_sta(:,iPG)*m_pesoPG_p(iPG); % Fuerzas internas debidas a la velocidad de filtracion del fluido

        velflu_new(:,iPG) = velflu_sta(:,iPG)-((mflu_new(:,iPG)/Dtime)*N4*(yy(1:2,1:4))')';
        velflu_total(:,iPG) = velflu_n1(:,iPG)-((mflu_new(:,iPG)/Dtime)*N4*(yy(1:2,1:4))')';

    elseif conshyp==16 %Modelo bifasico con solido elastico (saturado) con metodo de 
        %multiplicadores de Lagrange en desplazamientos. Caso Multi-escala
        m_sig_eps(:,:,iPG) = ct;   
        
        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*ct*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + B'*BiotM*N4*m_pesoPG_p(iPG); %Matriz de acople
        Qps = Qps + N4'*BiotM'*B*m_pesoPG_p(iPG); %Matriz de acople %Qps = Qps + N4*BiotM'*B*m_pesoPG_p(iPG);
        Sww = Sww + beta*(N4'*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); % Matriz de permeabilidad %Hww = Hww + DerivN*PermK*DerivN'*m_pesoPG_p(iPG);
        Llambda_u = Llambda_u + N8'*m_pesoPG_p(iPG); %Vector adicional para metodo de Mult. Lagrange por desplazamientos
        
        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); % Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velflu_sta(:,iPG)*m_pesoPG_p(iPG); % Fuerzas internas debidas a la velocidad de filtracion del fluido

        velflu_new(:,iPG) = velflu_sta(:,iPG)-((mflu_new(:,iPG)/Dtime)*N4*(yy(1:2,1:4))')';
        velflu_total(:,iPG) = velflu_n1(:,iPG)-((mflu_new(:,iPG)/Dtime)*N4*(yy(1:2,1:4))')';
        
    elseif conshyp==17  %Modelo bifasico con solido elastico (saturado) con metodo de 
        %multiplicadores de Lagrange en desplazamientos y poro presiones. Caso Multi-escala
        m_sig_eps(:,:,iPG) = ct;   
        
        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*ct*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + B'*BiotM*N4*m_pesoPG_p(iPG); %Matriz de acople
        Qps = Qps + N4'*BiotM'*B*m_pesoPG_p(iPG); %Matriz de acople %Qps = Qps + N4*BiotM'*B*m_pesoPG_p(iPG);
        Sww = Sww + beta*(N4'*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad
        Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); % Matriz de permeabilidad %Hww = Hww + DerivN*PermK*DerivN'*m_pesoPG_p(iPG);
        Llambda_u = Llambda_u + N8'*m_pesoPG_p(iPG); %Vector adicional para metodo de Mult. Lagrange por desplazamientos
        Llambda_p = Llambda_p + N4'*m_pesoPG_p(iPG); %Vector adicional para metodo de Mult. Lagrange por por opresiones
        
        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); % Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velflu_sta(:,iPG)*m_pesoPG_p(iPG); % Fuerzas internas debidas a la velocidad de filtracion del fluido

        velflu_new(:,iPG) = velflu_sta(:,iPG)-((mflu_new(:,iPG)/Dtime)*N4*(yy(1:2,1:4))')';
        velflu_total(:,iPG) = velflu_n1(:,iPG)-((mflu_new(:,iPG)/Dtime)*N4*(yy(1:2,1:4))')';
        
    elseif conshyp==50 %Modelo multiescala medio poroso saturado  
        m_sig_eps(:,:,iPG)= e_TanOp.m_sig_eps;
        m_chi_eps(:,:,iPG)= e_TanOp.m_chi_eps;
        m_V_eps(:,:,iPG)= e_TanOp.m_V_eps;
        m_sig_p(:,:,iPG)= e_TanOp.m_sig_p;
        m_chi_p(:,:,iPG)= e_TanOp.m_chi_p;
        m_V_p(:,:,iPG)= e_TanOp.m_V_p;
        m_sig_phi(:,:,iPG)= e_TanOp.m_sig_phi;
        m_chi_phi(:,:,iPG)= e_TanOp.m_chi_phi;
        m_V_phi(:,:,iPG)= e_TanOp.m_V_phi;

        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*m_sig_eps(:,:,iPG)*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + (B'*m_sig_phi(:,:,iPG)*DerivN+B'*m_sig_p(:,:,iPG)*N4)*m_pesoPG_p(iPG); %Matriz de acoplez
        Qps = Qps + (-DerivN'*m_V_eps(:,:,iPG)*B+N4'*m_chi_eps(:,:,iPG)*B)*m_pesoPG_p(iPG); %Matriz de acople
        Kpp = Kpp + (-DerivN'*m_V_phi(:,:,iPG)*DerivN-DerivN'*m_V_p(:,:,iPG)*N4+...
        N4'*m_chi_phi(:,:,iPG)*DerivN+N4'*m_chi_p(:,:,iPG)*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad + Matriz de permeabilidad!!!!?????       

        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); %Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); %Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velflu_new(:,iPG)*m_pesoPG_p(iPG); %Fuerzas internas debidas a la velocidad estatica de filtracion del fluido

    elseif conshyp==60 || conshyp==61 %Modelo multiescala medio poroso saturado - Resolucion Analitica  
        m_sig_eps(:,:,iPG)= e_TanOp.m_sig_eps;
        m_chi_eps(:,:,iPG)= e_TanOp.m_chi_eps;
        m_V_eps(:,:,iPG)= e_TanOp.m_V_eps;
        m_sig_p(:,:,iPG)= e_TanOp.m_sig_p;
        m_chi_p(:,:,iPG)= e_TanOp.m_chi_p;
        m_V_p(:,:,iPG)= e_TanOp.m_V_p;
        m_sig_phi(:,:,iPG)= e_TanOp.m_sig_phi;
        m_chi_phi(:,:,iPG)= e_TanOp.m_chi_phi;
        m_V_phi(:,:,iPG)= e_TanOp.m_V_phi;

        % Calculo de las sumbatrices del la matriz de rigidez acoplada
        Kss = Kss + B'*m_sig_eps(:,:,iPG)*B*m_pesoPG_d(iPG); %Matriz de rigidez
        Qsp = Qsp + (B'*m_sig_phi(:,:,iPG)*DerivN+B'*m_sig_p(:,:,iPG)*N4)*m_pesoPG_p(iPG); %Matriz de acoplez
        Qps = Qps + (-DerivN'*m_V_eps(:,:,iPG)*B+N4'*m_chi_eps(:,:,iPG)*B)*m_pesoPG_p(iPG); %Matriz de acople
        Kpp = Kpp + (-DerivN'*m_V_phi(:,:,iPG)*DerivN-DerivN'*m_V_p(:,:,iPG)*N4+...
        N4'*m_chi_phi(:,:,iPG)*DerivN+N4'*m_chi_p(:,:,iPG)*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad + Matriz de permeabilidad!!!!?????       

        % Calculo de fuerzas internas (P/residuo)
        fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); %Fuerzas internas debidas a las tensiones totales
        fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); %Fuerzas internas debidas al contenido de masa del fluido
        fflux_V = fflux_V + DerivN'*velflu_new(:,iPG)*m_pesoPG_p(iPG); %Fuerzas internas debidas a la velocidad estatica de filtracion del fluido
    end %if(conshyp)   
end %for(ipg)

if conshyp==14 %Caso mono-escala
    %Segun sea el argumento Full Order Expande o Zero Order Expended
    Beta_Sww = beta_factor2*Sww;
    Kpp = Beta_Sww + theta*Dtime*Hww;%Matriz de 'rigidez' al flujo

    fflux  = fflux - fflux_m + Dtime*fflux_V; %Fuerzas internas asociadas al flujo

    % Vector de fuerzas elemental acoplado
    Fe = [fstreT ; fflux ; zeros(4,1)];
    
    %Segun sea el argumento Full Order Expande o Zero Order Expended
    Beta_Qsp = beta_factor2*Qsp;
    
    % Matriz de rigidez elemental acoplada
    Ke = [Kss -Beta_Qsp zeros(16,4);
         -Qps -Kpp zeros(4,4) ;
          zeros(4,20) eye(4,4)];
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand(Ke,Fe,dofpe) ;
elseif conshyp==15
    %Segun sea el argumento Full Order Expande o Zero Order Expended
    Beta_Sww = beta_factor2*Sww;
    Kpp = Beta_Sww + Dtime*theta*Hww;%Matriz de 'rigidez' al flujo
        
    %Fuerzas internas asociadas al flujo con multiplicadores de Lagrange en
    %poro presiones
    fflux  = fflux - fflux_m + Dtime*fflux_V + Llambda_p*dulambda_p;
    
    %Fuerzas internas asociadas al multiplicador de Lagrange en poro
    %presiones
    fMul_Lag_p = fMul_Lag_p + Llambda_p'*dup;
    
    % Vector de fuerzas elemental acoplado
    Fe = [fstreT ; fflux ; zeros(4,1); fMul_Lag_p];
   
    %Segun sea el argumento Full Order Expande o Zero Order Expended
    Beta_Qsp = beta_factor2*Qsp;
    
    %Factor para la diagonal de la matriz expandida 
    %por incognitas adicionales debidas a los mult. Lagrange
    fd_p = 0.00000000000001/e_VG.nElem;

    
    % Matriz de rigidez elemental acoplada
    Ke = [Kss -Beta_Qsp zeros(16,5);
         -Qps -Kpp zeros(4,4) Llambda_p;
          zeros(4,20) eye(4,4) zeros(4,1);
          zeros(1,16) Llambda_p' zeros(1,4) fd_p];
      
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand(Ke,Fe,dofpe+1) ;
        
elseif conshyp==16
    %Segun sea el argumento Full Order Expande o Zero Order Expended
    Beta_Sww = beta_factor2*Sww;
    Kpp = Beta_Sww + Dtime*theta*Hww;%Matriz de 'rigidez' al flujo
        
   
    fMul_Lag_u = fMul_Lag_u - Llambda_u'*dud; 
%     fMul_Lag_u = fMul_Lag_u - Llambda_u'*Dtime*(theta*dud+ud_old); 

    %Fuerzas internas asociadas al flujo con multiplicadores de Lagrange en
    %poro presiones
    fflux  = fflux - fflux_m + Dtime*fflux_V;
    
    % Vector de fuerzas elemental acoplado
    Fe = [fstreT ; fflux ; zeros(4,1); fMul_Lag_u];
   
    %Segun sea el argumento Full Order Expande o Zero Order Expended
    Beta_Qsp = beta_factor2*Qsp;
    
    %Factor para la diagonal de la matriz expandida 
    %por incognitas adicionales debidas a los mult. Lagrange
    fd_d = 0.00000000000001/e_VG.nElem; 
     
%     Llambda_u = Dtime*theta*Llambda_u;

    
    % Matriz de rigidez elemental acoplada
    Ke = [Kss -Beta_Qsp zeros(16,4) -Llambda_u;
         -Qps -Kpp zeros(4,6) ;
          zeros(4,20) eye(4,4) zeros(4,2);
          -Llambda_u' zeros(2,8) fd_d*eye(2,2)];
      
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand(Ke,Fe,dofpe+2) ;
    
elseif conshyp==17
    %Segun sea el argumento Full Order Expande o Zero Order Expended
    Beta_Sww = beta_factor2*Sww;
    Kpp = Beta_Sww + Dtime*theta*Hww;%Matriz de 'rigidez' al flujo
        
    %Fuerzas internas asociadas a las tensiones con con multiplicadores de Lagrange en
    %desplazamientos
    fstreT = fstreT - Llambda_u*dulambda_u;
%     fstreT = fstreT - Llambda_u*Dtime*(theta*dulambda_u+ulambda_u);
    
    fMul_Lag_u = fMul_Lag_u - Llambda_u'*dud; 
%     fMul_Lag_u = fMul_Lag_u - Llambda_u'*Dtime*(theta*dud+ud_old); 

    %Fuerzas internas asociadas al flujo con multiplicadores de Lagrange en
    %poro presiones
    fflux  = fflux - fflux_m + Dtime*fflux_V + Llambda_p*dulambda_p;
%     fflux  = fflux - fflux_m + Dtime*fflux_V + Llambda_p*Dtime*(theta*dulambda_p+ulambda_p);
    
    %Fuerzas internas asociadas al multiplicador de Lagrange en poro
    %presiones
    fMul_Lag_p = fMul_Lag_p + Llambda_p'*dup;
%     fMul_Lag_p = fMul_Lag_p + Llambda_p'*Dtime*(theta*dup+up_old); 
    
    % Vector de fuerzas elemental acoplado
    Fe = [fstreT ; fflux ; zeros(4,1); fMul_Lag_u; fMul_Lag_p];
   
    %Segun sea el argumento Full Order Expande o Zero Order Expended
    Beta_Qsp = beta_factor2*Qsp;
    
    %Factor para la diagonal de la matriz expandida 
    %por incognitas adicionales debidas a los mult. Lagrange
    fd_d = 0.00000000000001/e_VG.nElem; 
    fd_p = 0.00000000000001/e_VG.nElem;
     
%     Llambda_u = Dtime*theta*Llambda_u;
%     Llambda_p = Dtime*theta*Llambda_p;
    
    % Matriz de rigidez elemental acoplada
    Ke = [Kss -Beta_Qsp zeros(16,4) -Llambda_u zeros(16,1);
         -Qps -Kpp zeros(4,6) Llambda_p;
          zeros(4,20) eye(4,4) zeros(4,3);
          -Llambda_u' zeros(2,8) fd_d*eye(2,2) zeros(2,1);
          zeros(1,16) Llambda_p' zeros(1,6) fd_p];
      
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand(Ke,Fe,dofpe+3) ;
    
elseif conshyp==50 || conshyp==60 || conshyp==61
    fflux  = fflux - fflux_m + Dtime*fflux_V; %Fuerzas internas asociadas al flujo    
%     fflux  = fflux - fflux_m + Dtime*fflux_V - fflux_mY; %Fuerzas internas asociadas al flujo    
    
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
%         Ke = [Kss  Qsp zeros(16,4);
%          -Qps -Kpp zeros(4,4) ;
%          zeros(4,20) eye(4,4)];
    Ke = [Kss*Dtime  Qsp*Dtime zeros(16,4);
         -Qps*Dtime -Kpp*Dtime zeros(4,4) ;
         zeros(4,20) eye(4,4)];
     
    % Expande y ordena la matriz de rigidez elemental
    [kt,fint] = f_Expand(Ke,Fe,dofpe) ;
end %if(conshyp) 

%############################################################################################
%SALIDA DE INFORMACION
% Se ordena las matrices como vectores columnas

%VARIABLES PRIMITIVAS
%Delta de deformaciones
eps_new    = eps_new(:); 
%Delta del gradiente de poro presiones
phi_new    = phi_new(:); 
%Delta de poro presiones
p_new    = p_new(:); 

%VARIABLES PRIMITIVAS FLUCTUANTES
%Delta de micro-deformaciones fluctuantes o delta de macro-deformaciones
eps_fluct  = eps_fluct(:); 
%Delta de gradiente de micro-poro presiones fluctuantes o delta de gradiente de macro-poro presiones
phi_fluct  = phi_fluct(:); 
%Delta de micro-poro presiones fluctuantes o delta de macro-poro presiones
p_fluct  = p_fluct(:); 

%VARIABLES DUALES
%Delta de tensiones efectivas
sigmaE_new  = sigmaE_new(:); 
%Delta de tensiones totales
sigmaT_new  = sigmaT_new(:);
%Delta de tasa del contenido de masa del fluido
mflu_new  = mflu_new(:);
%Componente estatica de la velocidad de filtracion  al tiempo "n+theta"
velflu_sta = velflu_sta(:);
% %Delta de la componente dinamica de la velocidad de filtracion
% velflu_dyn = velflu_dyn(:);

%Velocidad de filtracion  al tiempo "n+theta"
velflu_new = velflu_new(:);

%Velocidad de filtracion  al tiempo "n+1"
velflu_total = velflu_total(:);

%Variables historicas de cada modelo constitutivo
hvar_new   = hvar_new(:);
aux_var    = aux_var(:);
