function [kt,fint,sigma_new,eps_new,eps_fluct,hvar_new,aux_var,m_TensorTang] = ...
   f_MatElem_bbar_q1(u,hvar_old,aux_var,e_DatElemSet,e_DatMatSet,...
   m_Be,m_DetJe,DefMacro,sigma_old,e_VG)

% Variables globales
ntens = e_VG.ntens;
dofpe = e_DatElemSet.dofpe;
nPG = e_DatElemSet.npg; 
wg = e_DatElemSet.wg;

sihvarpg = e_DatMatSet.sihvarpg;
siavarpg = e_DatMatSet.siavarpg;
conshyp  = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;
      
% Inicializaciones
kt           = zeros(dofpe,dofpe);
fint         = zeros(dofpe,1);
sigma_new    = zeros(ntens,nPG);
eps_new      = zeros(ntens,nPG);
eps_fluct    = zeros(ntens,nPG);
%hvar_new     = zeros(sihvarpg,nPG);
hvar_new = f_DefinicionhVar(conshyp,sihvarpg,nPG);
%aux_var      = zeros(siavarpg,nPG);
if esImplex
   m_TensorTang = zeros(ntens,ntens,2*nPG);
else
   m_TensorTang = zeros(ntens,ntens,nPG);
end
hvar_old     = reshape(hvar_old,sihvarpg,[]);
%hvar_new = reshape(hvar_new,sihvarpg,[]);
aux_var = reshape(aux_var,siavarpg,nPG);
%sigma_old    = reshape(sigma_old,ntens,nPG);
%DefMacro     = reshape(DefMacro,ntens,[]);
%m_pesoPG     = m_DetJe.*wg*thickness;
%En el peso de gauss ya viene multiplicado el espesor.
m_pesoPG = m_DetJe.*wg;
%volElem = sum(m_pesoPG);
 
for iPG = 1:nPG
   
   e_VG.iPG = iPG;

   B = m_Be(:,:,iPG);

   eps_fluct(:,iPG) = B*u;

   % Deformacion aplicada a cada punto de Gauss
%     if e_VG.esME
        eps_new(:,iPG) = DefMacro + eps_fluct(:,iPG); %JLM
%     else    
%        eps_new(:,iPG) = DefMacro(:,iPG)+eps_fluct(:,iPG); % AA: Bloque� linea 50
%     end
            
   % Modelo constitutivo
   switch conshyp
      case 1   %Elasticidad
         [ct,sigma_new(:,iPG)] = rmap_elast(eps_new(:,iPG),e_DatMatSet.ce);
      case 2   %ELASTO - PLASTICIDAD J2: HARDENING-SOFTENING ISOTROPICO
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 4   %Da�o is�tropico
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG); 
      case 9
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2bl(...
            ntens,eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 10  %Da�o isotr�pico regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_reg(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 11  %Da�o isotr�pico solo tracci�n regularizado
          if e_VG.elast                           
              [ct,sigma_new(:,iPG),~,hvar_new(:,iPG)] ...
                          = rmap_damage_elastic(eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet);
          else
              Eprop_qSD.ksd = e_DatMatSet.m_hReg(e_VG.iElemSet);
              [ct,sigma_new(:,iPG),~,hvar_new(:,iPG)] = rmap_damage_IMPLEX (eps_new(:,iPG),...
                  hvar_old(:,iPG),e_DatMatSet,Eprop_qSD);
          end
      case 12  %Da�o isotr�pico solo tracci�n regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] =...
             RMapDanoSTraccReg (...
               eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      case 13 %Da�o isotr�pico solo tracci�n regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = ...
            RMapDanoRankine(...
            eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      case 50  %Modelo multiescala cl�sico          
%          %fprintf('**-- Inicio del return mapping del modelo multiescala\n')
%          cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
%          listeners = cmdWinDoc.getDocumentListeners;
%          jFxCommandArea = listeners(3);
%          %colorTextArea = get(jTextArea,'Background');
%          set(jFxCommandArea,'Background','red');
         [ct,sigma_new(:,iPG),hvar_new(:,iPG)] = f_RMap_ME(eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
%          set(jFxCommandArea,'Background','yellow');
%          %fprintf('**-- Fin del return mapping del modelo multiescala\n')
      case 52  %Elasticidad usando el tensor el�stico obtenido de una homogenizaci�n de una microcelda.
         [ct,sigma_new(:,iPG)] = rmap_elast(eps_new(:,iPG),e_DatMatSet.ce);
      otherwise
         error('Matrices Elementales BBar_q1: Modelo constitutivo no definido.')
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
   
   % C�lculo de fint
   %fint = fint+B'*(sigma_new(:,iPG)-sigma_old(:,iPG))*m_pesoPG(iPG);
   fint = fint+B'*sigma_new(:,iPG)*m_pesoPG(iPG);  
   
   % C�lculo de Kep
   kt = kt+B'*ct*B*m_pesoPG(iPG);
     
      
   % C�lculo de la energia interna
   switch conshyp
       case 1 % Elastic material
           IntEnergyNew(iPG) = 0;
       case 2 % Plastic material W = 0.5*(e_e*C*e_e)+0.5*K*chi^2
           IntEnergyNew(iPG) = 0.5*((eps_new(:,iPG)-hvar_new(1:4,iPG))'*e_DatMatSet.ce*(eps_new(:,iPG)-hvar_new(1:4,iPG)))...
               +0.5*hvar_new(6,iPG)*(hvar_new(5,iPG))^2; % PRUEBA!!
%                +0.5*hvar_new(12,iPG)*(hvar_new(5,iPG))^2; % PRUEBA!!           
       case 11  % Continuum damage model
           IntEnergyNew(iPG) = 0.5*((1-aux_var(1,iPG))*eps_new(:,iPG)'*...
               e_DatMatSet.ce*eps_new(:,iPG));  % energia interna
       case 50 % Multiscale classic model
           %AA: desbloque la linea de abajo y bloque la que NO ENTIENDO   
           IntEnergyNew(iPG) = 0.5*(sigma_new(:,iPG)'*eps_new(:,iPG));  % energia interna
           
           % Flag to split the snapshots (elastic & inelastic part)
           %AA: NO ENTIENDO. Quiere tomar el valor de la deformaci�n
           % pl�stica como valor de la energia interna???
%            IntEnergyNew(iPG) = max(max(hvar_new(1).e_VarEst(1).hvar(13:e_DatMatSet.e_DatSet(1).e_DatMat.sihvarpg:end,:)));

       otherwise
           IntEnergyNew(iPG) = 0.5*((1-aux_var(1,iPG))*eps_new(:,iPG)'*...
               e_DatMatSet.ce*eps_new(:,iPG));  % energia interna
   end
   
   
   
   
   
   
   
   
   
   
end

% Se ordena las matrices como vectores columnas
sigma_new  = sigma_new(:);
hvar_new   = hvar_new(:);
aux_var    = aux_var(:);
eps_new    = eps_new(:);
eps_fluct  = eps_fluct(:);
   
end