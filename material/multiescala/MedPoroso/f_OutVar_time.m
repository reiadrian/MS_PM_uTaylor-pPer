%Funcion que determina las variables en el tiempo "n+1" o "n+theta"
%segun la necesidad
function e_VarEst_new = f_OutVar_time(e_VarEst_new,e_VarEst_old,...
            eps_new,phi_new,porpr_new,eps_fluct,phi_fluct,p_fluct,sigmaE_new,...
            sigmaT_new,mflu_new,velflu_sta,velflu_total,velflu_new,e_VG,flag)

if flag==0
%############################################################################################
%VARIABLES EN EL PASO DE TIEMPO "n+1"
%############################################################################################
theta = e_VG.theta;
Dtime = e_VG.Dtime;
%VARIABLES PRIMITIVAS
%AL PASO DE TIEMPO "n+1"
%Deformaciones totales al paso de tiempo "n+1"
e_VarEst_new.eps = e_VarEst_old.eps + eps_new;
%Gradiente de poro presiones totales al paso de tiempo "n+1"
e_VarEst_new.phi =  e_VarEst_old.phi + phi_new;              
%Poro presiones  totales al paso de tiempo "n+1"
e_VarEst_new.porpr = e_VarEst_old.porpr + porpr_new;

%VARIABLES PRIMITIVAS FLUCTUANTES
%Deformaciones fluctuantes al paso de tiempo "n+1"
e_VarEst_new.eps_fluct = e_VarEst_old.eps_fluct + eps_fluct; 
%Gradiente de poro presiones fluctuantes al paso de tiempo "n+1"
e_VarEst_new.phi_fluct = e_VarEst_old.phi_fluct + phi_fluct;
%Poro presiones fluctuantes al paso de tiempo "n+1"
e_VarEst_new.p_fluct = e_VarEst_old.p_fluct + p_fluct;

%VARIABLES DUALES
%Tensiones efectivas al paso de tiempo "n+1"
e_VarEst_new.sigmaE = e_VarEst_old.sigmaE + sigmaE_new;
%Tensiones totales al paso de tiempo "n+1"
e_VarEst_new.sigmaT = e_VarEst_old.sigmaT + sigmaT_new;
%Tasa de contenido de masa del fluido al paso de tiempo "n+1"
e_VarEst_new.mflu = e_VarEst_old.mflu + mflu_new;        

if e_VG.conshyp == 50 || e_VG.conshyp == 60 || e_VG.conshyp == 61 || e_VG.conshyp == 62
%     %Componente estatica de la velocidad de filtracion  al tiempo "n+1" 
%     e_VarEst_new.velflu_sta = velflu_sta;
% %     e_VarEst_new.velflu_sta = (velflu_sta/theta)-((1-theta)*e_VarEst_old.velflu_sta/theta);
%     %Componente dinamica de la velocidad de filtracion  al tiempo "n+1" 
%     e_VarEst_new.velflu_dyn = velflu_dyn; 
% %     e_VarEst_new.velflu_dyn = velflu_dyn  + e_VarEst_old.velflu_dyn; 
%     %Velocidad de filtracion homogeneizada al paso de tiempo "n+theta" 
    e_VarEst_new.velflu = velflu_new;
%     %Velocidad de filtracion homogeneizada al paso de tiempo "n+1" 
    e_VarEst_new.velflu_total = velflu_total;
% %     velflu_i= e_VarEst_new.velflu_sta - e_VarEst_new.velflu_dyn; 
% %     e_VarEst_new.velflu = e_VarEst_new.velflu_sta - e_VarEst_new.velflu_dyn ;
%     e_VarEst_new.velflu = (velflu_total/theta)-((1-theta)*e_VarEst_old.velflu/theta);
elseif isfield(e_VG,'conshypMacro') && e_VG.conshyp ~= 50 && e_VG.conshyp ~= 60 && e_VG.conshyp ~= 61
    %Velocidad de filtracion al paso de tiempo "n+1" 
    e_VarEst_new.velflu = velflu_sta;
    e_VarEst_new.velflu_total = (velflu_sta/theta)-((1-theta)*e_VarEst_old.velflu/theta);
elseif ~isfield(e_VG,'conshypMacro') && e_VG.conshyp==14
    e_VarEst_new.velflu = velflu_new;
    e_VarEst_new.velflu_total = velflu_total;    
end

elseif flag==1
    %############################################################################################
    %VARIABLES MACRO-ESCALA EN EL PASO DE TIEMPO "n+theta"
    %############################################################################################
    theta = e_VG.theta;
     %VARIABLES PRIMITIVAS
    %AL PASO DE TIEMPO "n+1"
    %Deformaciones totales al paso de tiempo "n+theta"
    e_VarEst_new.eps = e_VarEst_old.eps + theta*eps_new;
    %Gradiente de poro presiones totales al paso de tiempo "n+theta"
    e_VarEst_new.phi = e_VarEst_old.phi + theta*phi_new;
    %Poro presiones  totales al paso de tiempo "n+theta"
    e_VarEst_new.porpr = e_VarEst_old.porpr + theta*porpr_new;

    %VARIABLES PRIMITIVAS FLUCTUANTES
    %Deformaciones fluctuantes al paso de tiempo "n+theta"
    e_VarEst_new.eps_fluct = e_VarEst_old.eps_fluct + theta*eps_fluct; 
    %Gradiente de poro presiones fluctuantes al paso de tiempo "n+theta"
    e_VarEst_new.phi_fluct = e_VarEst_old.phi_fluct + theta*phi_fluct;
    %Poro presiones fluctuantes al paso de tiempo "n+theta"
    e_VarEst_new.p_fluct = e_VarEst_old.p_fluct + theta*p_fluct;

    %VARIABLES DUALES
    %Tensiones efectivas al paso de tiempo "n+theta"
    e_VarEst_new.sigmaE = e_VarEst_old.sigmaE + theta*sigmaE_new;
    %Tensiones totales al paso de tiempo "n+theta"
    e_VarEst_new.sigmaT = e_VarEst_old.sigmaT + theta*sigmaT_new;
    %Tasa de contenido de masa del fluido al paso de tiempo "n+theta"
    e_VarEst_new.mflu = e_VarEst_old.mflu + theta*mflu_new;        

    if e_VG.conshyp == 50
        %Componente estatica de la velocidad de filtracion  al tiempo "n+theta" 
        e_VarEst_new.velflu_sta = velflu_sta;
        %Componente dinamica de la velocidad de filtracion  al tiempo "n+theta" 
        e_VarEst_new.velflu_dyn =  e_VarEst_old.velflu_dyn + theta*velflu_dyn; 
        %Velocidad de filtracion homogeneizada al paso de tiempo "n+theta" 
        e_VarEst_new.velflu = e_VarEst_new.velflu_sta - e_VarEst_new.velflu_dyn;
    else
        %Velocidad de filtracion al paso de tiempo "n+theta" 
        e_VarEst_new.velflu = velflu_sta;
    end
elseif flag==2
        %############################################################################################
        %VARIABLES MICRO-ESCALA (SOLAMENTE NECESARIOS PARA LA MICRO)
        %############################################################################################
        %VARIABLES PRIMITIVAS
        %Delta de deformaciones
        e_VarEst_new.eps = eps_new;
        %Delta del gradiente de poro presiones
        e_VarEst_new.phi = phi_new;               
        %Delta de poro presiones
        e_VarEst_new.porpr = porpr_new;

        %VARIABLES PRIMITIVAS FLUCTUANTES
        %Delta de deformaciones fluctuantes
        e_VarEst_new.eps_fluct = eps_fluct; 
        %Delta de gradiente de poro presiones fluctuantes
        e_VarEst_new.phi_fluct = phi_fluct;        
        %Delta de poro presiones fluctuantes
        e_VarEst_new.p_fluct = p_fluct;

        %VARIABLES DUALES
        %############################################################################################
        %PERMITEN CALCULAR VARIABLES DUALES MACRO-ESCALA HOMOGENEIZADAS (NECESARIO 
        %LLEVAR EN CUENTA TANTO PARA LA MICRO-ESCALA COMO LA MACRO-ESCALA)
        %############################################################################################
        %Delta de tensiones efectivas
        e_VarEst_new.sigmaE = sigmaE_new;
        %Delta de tensiones totales
        e_VarEst_new.sigmaT = sigmaT_new;
        %Deltra de tasa de contenido de masa del fluido
        e_VarEst_new.mflu = mflu_new;        
        %Componente estatica de la velocidad de filtracion  al tiempo "n+theta" - UTIL
        %EN LA MICRO-ESCALA A FINES DE HOMOGENEIZAR
        e_VarEst_new.velflu_sta = velflu_sta; 
        %Delta de la componente dinamica de la velocidad de filtracion - UTIL
        %EN LA MICRO-ESCALA A FINES DE HOMOGENEIZAR
        %Velocidad de filtracion al paso de tiempo "n+1" 
        e_VarEst_new.velflu_total = velflu_total; 
         %Velocidad de filtracion al paso de tiempo "n+theta" 
        e_VarEst_new.velflu = velflu_new;
    
end
