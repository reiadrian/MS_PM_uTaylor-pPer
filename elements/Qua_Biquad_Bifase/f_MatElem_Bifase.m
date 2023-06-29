function [kt,fint,sigmaE_new,sigmaT_new,eps_new,eps_fluct,hvar_new,aux_var,m_TensorTang] = ...
   f_MatElem_Bifase(u,eps_old,hvar_old,aux_var,e_DatElemSet,e_DatMatSet,...
   m_Be_d,m_DetJe_d,m_Dercae_p,m_DetJe_p,m_FF_p,DefMacro,sigmaE_old,sigmaT_old,e_VG)

%AA: Cree funci�n

% Variable globales
ntens    = e_VG.ntens;
theta    = e_VG.theta;
DTime    = e_VG.Dtime;
%ndn    = e_VG.ndn;
%ndime    = e_VG.ndime;
dofpe = e_DatElemSet.dofpe;
%npe = e_DatElemSet.npe;
dofpe_d = e_DatElemSet.dofpe_d; %AA
dofpe_p = e_DatElemSet.dofpe_p; %AA
nPG = e_DatElemSet.npg; 
wg = e_DatElemSet.wg;
pos_d = e_DatElemSet.pos_d;
pos_p = e_DatElemSet.pos_p;

sihvarpg = e_DatMatSet.sihvarpg;
siavarpg = e_DatMatSet.siavarpg;
conshyp  = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;

% Inicializaciones
Kss           = zeros(dofpe_d,dofpe_d);
Qsp           = zeros(dofpe_d,dofpe_p);
Qps           = zeros(dofpe_p,dofpe_d);
Sww           = zeros(dofpe_p,dofpe_p);
Hww           = zeros(dofpe_p,dofpe_p);
fstreT        = zeros(dofpe_d,1); %Fuerzas internas debidas a las tensiones totales

fflux_m         = zeros(dofpe_p,1); %Fuerzas internas asociadas al flujo
fflux_V         = zeros(dofpe_p,1); %Fuerzas internas asociadas al flujo

fflux         = zeros(dofpe_p,1); %Fuerzas internas asociadas al flujo

fint          = zeros(dofpe,1); % VER QUE HACER CON LAS FUERZAS!!!!
sigmaE_new     = zeros(ntens,nPG); %AA
sigmaT_new     = zeros(ntens,nPG); %AA
eps_new       = zeros(ntens,nPG);

phi_new       = zeros(2,nPG);
velflu_new    = zeros(2,nPG);  % Corresponde a PG del elemento en la MICRO
mflu_new    = zeros(1,nPG); % Corresponde a PG del elemento en la MICRO


eps_fluct     = zeros(ntens,nPG);
%hvar_new     = zeros(sihvarpg,nPG);
hvar_new      = f_DefinicionhVar(conshyp,sihvarpg,nPG); %AA: add case 14
%aux_var      = zeros(siavarpg,nPG);

if esImplex
   m_TensorTang = zeros(ntens,ntens,2*nPG);
else
   m_TensorTang = zeros(ntens,ntens,nPG);
end

hvar_old = reshape(hvar_old,sihvarpg,[]);
%hvar_new = reshape(hvar_new,sihvarpg,[]);
aux_var = reshape(aux_var,siavarpg,nPG);
%sigma_old = reshape(sigma_old,ntens,nPG);
%eps_old = reshape(eps_old,ntens,[]);
%DefMacro = reshape(DefMacro,ntens,[]);

%En el peso de gauss ya viene multiplicado el espesor.
m_pesoPG_d = m_DetJe_d.*wg;
m_pesoPG_p = m_DetJe_p.*wg;

ud = u(pos_d);
up = u(pos_p);

for iPG = 1:nPG
   
   e_VG.iPG = iPG;

   B = m_Be_d(:,:,iPG);
   
   DerivN = m_Dercae_p(:,:,iPG);
   N4 = m_FF_p(:,:,iPG);
   
   eps_fluct(:,iPG) = B*ud;

   % Deformacion aplicada a cada punto de Gauss
   %eps_new(:,iPG) = DefMacro(:,iPG)+eps_fluct(:,iPG);
   eps_new(:,iPG) = DefMacro + eps_fluct(:,iPG);
   pNodo_new = up;
   phi_new(:,iPG) = DerivN*up;

   % Modelo constitutivo
   switch conshyp
       case 14          
         [ct,BiotM,beta,PermK,sigmaE_new(:,iPG),sigmaT_new(:,iPG),mflu_new(1,iPG),velflu_new(:,iPG)] = ...
                f_Rmap_Bif_ElastME(eps_new(:,iPG),phi_new(:,iPG),pNodo_new,N4,e_DatMatSet);  %AA: cree funcion
%          [ct,BiotM,beta,PermK,sigmaE_new(:,iPG),sigmaT_new(:,iPG)] = ...
%              f_Rmap_Bif_Elast(eps_new(:,iPG),e_DatMatSet,up,N4);  %AA: cree funci�n
%        case {50,55}  %Modelo multiescala cl�sico
%          %fprintf('**-- Inicio de return mapping del modelo multiescala\n')
%          %cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
%          %listeners = cmdWinDoc.getDocumentListeners;
%          %jFxCommandArea = listeners(3);
%          %set(jFxCommandArea,'Background','red');
%          [ct,sigma_new(:,iPG),hvar_new(:,iPG)] = f_RMap_ME(eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
%          %set(jFxCommandArea,'Background','yellow');
%          %fprintf('**-- Fin de return mapping del modelo multiescala\n')
%        case 52  %Elasticidad usando el tensor el�stico obtenido de una homogenizaci�n de una microcelda.
%          [ct,sigma_new(:,iPG)] = rmap_elast(eps_new(:,iPG),e_DatMatSet.ce);
%        case 110 %Large deformations J2 Plasticity
%          [ct,sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG),strsg,dmatx] = ...
%             f_RMapPlastJ2LD(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      otherwise
         error('Matrices Elementales bifasicas: Modelo constitutivo no definido.')
   end

   %Se almacena para los tensor tangente constitutivos para realizar homogeneizaci�n y an�lisis de
   %bifurcaci�n
   if esImplex
      % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
      %tangentes impl�citas para el an�lisis de bifurcaci�n.
      %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando est� seleccionado
      %el implex.
      m_TensorTang(:,:,iPG) = ct.Implex;
      %Se almacena los tensores tangentes constitutivo impl�citos para an�lisis de bifurcaci�n como si fuera
      %PG adicionales, tantos como nPG. Se almacena en los �ndices (:,:,nPG+1:2*nPG).
      %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
      %tercera dimensi�n de m_TensorTang (size(m_TensorTang,3)).
      m_TensorTang(:,:,iPG+nPG) = ct.Impli;
      %En los c�lculos para el ensamblaje se utiliza el implex.
      ct = ct.Implex;
   else
      if isfield(ct,'Impli')  
         ct = ct.Impli;
      end
      m_TensorTang(:,:,iPG) = ct;
   end
   
  
   % C�lculo de las sumbatrices del la matriz de rigidez acomplada
   Kss = Kss + B'*ct*B*m_pesoPG_d(iPG); %Matriz de rigidez
   Qsp = Qsp + B'*BiotM*N4*m_pesoPG_p(iPG); %Matriz de acople
   Qps = Qps + N4'*BiotM'*B*m_pesoPG_p(iPG); %Matriz de acople %Qps = Qps + N4*BiotM'*B*m_pesoPG_p(iPG);
   Sww = Sww + beta*(N4'*N4)*m_pesoPG_p(iPG); %Matriz de compresibilidad
   Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); % Matriz de permeabilidad %Hww = Hww + DerivN*PermK*DerivN'*m_pesoPG_p(iPG);
   
%    % C�lculo de fuerzas internas (P/residuo)
%    fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); %Fuerzas internas debidas a las tensiones totales
% Calculo de fuerzas internas (P/residuo)
fstreT = fstreT + B'*sigmaT_new(:,iPG)*m_pesoPG_d(iPG); % Fuerzas internas debidas a las tensiones totales
%mflu_new = BiotM'*eps_new + beta*N4*p_new; % Contenido de masa del fluido
fflux_m = fflux_m + N4'*mflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas al contenido de masa del fluido"
%velflu_new = -PermK*DerivN*phi_new; % Velocidad de filtracion del fluido
fflux_V = fflux_V + DerivN'*velflu_new(:,iPG)*m_pesoPG_p(iPG); % "Fuerzas internas debidas a la velocidad de filtracion del fluido"
       
end %for(ipg)

   Kpp = Sww + (theta*DTime*Hww);
   fflux  = fflux - fflux_m - (theta*DTime)*fflux_V; %Fuerzas internas asociadas al flujo
%    fflux  = fflux - Qps*ud - Kpp*up; %Fuerzas internas asociadas al flujo
   %fflux  = fflux - N4'*BiotM'*B*ud - beta*N4'*N4*up + theta*DTime*DerivN'*PermK*DerivN*up; %Fuerzas internas asociadas al flujo
   
   % Matriz de rigidez elemental acoplada
     Ke = [Kss -Qsp zeros(16,4);
        -Qps -Kpp zeros(4,4) ;
        zeros(4,20) eye(4,4)];
    
    % Vector de fuerzas elemental acoplado
      Fe = [fstreT ; fflux ; zeros(4,1)];
    
    % Expande y ordena la matriz de rigidez elemental
      [kt,fint] = f_Expand(Ke,Fe,dofpe) ; 
       
% Se ordena las matrices como vectores columnas
sigmaE_new  = sigmaE_new(:); %AA (por PG)
sigmaT_new  = sigmaT_new(:); %AA (por PG)
hvar_new   = hvar_new(:);
aux_var    = aux_var(:);
eps_new    = eps_new(:);
eps_fluct  = eps_fluct(:);
   
