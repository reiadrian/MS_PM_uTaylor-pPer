function [kt,fint,sigma_new,eps_new,eps_fluct,hvar_new,aux_var,m_TensorTang,...
    IntEnergyNew] = f_MatElem_FBar_q1(u,hvar_old,aux_var,m_VarAuxElem,...
    e_DatElemSet,e_DatMatSet,m_Be,m_DetJe,DefMacro,eps_old,sigma_old,e_VG)

% Variable globales
%conshyp = e_VG.conshyp;
%npe = e_VG.npe;
%dofpe = e_VG.dofpe;
ntens = e_VG.ntens;
%struhyp = e_VG.struhyp;
%
dofpe = e_DatElemSet.dofpe;
nPG = e_DatElemSet.npg; 
wg = e_DatElemSet.wg;
sihvarpg = e_DatMatSet.sihvarpg;
siavarpg = e_DatMatSet.siavarpg;
conshyp  = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;

% Inicializaciones
kt = zeros(dofpe,dofpe);
fint = zeros(dofpe,1);
sigma_new = zeros(ntens,nPG);
eps_new = zeros(ntens,nPG);
eps_fluct = zeros(ntens,nPG);
IntEnergyNew  = zeros(nPG,1) ;  % el 4 es para llevar 4 energias, Psi_e, Psi_e_vol, Psi_e_dev, Psi_p
hvar_new = f_DefinicionhVar(conshyp,sihvarpg,nPG);
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
%En el peso de gauss ya viene multiplicado el espesor.
m_pesoPG = m_DetJe.*wg;

% Punto de gauss central.
%Determinación del determinante de gradiente de deformación del punto central
B = reshape(m_VarAuxElem,ntens,[]);
%En DefMacro viene la identidad en el caso de los problemas monoescala.
m_F0 = DefMacro+B*u;
%Matriz 2D de F0
m_F2D = [m_F0(1),m_F0(4);m_F0(5),m_F0(2)];
Fzz = m_F0(3);
%Se utiliza diferentes formas de calcular el FBar para deformación plana.
tipoFBar = 1;
switch tipoFBar
   case 1  %Plane deformation
      %Se utiliza lo propuesto por Souza, de no modificar la componente Fzz del tensor.
      J0 = det(m_F2D);
      pot = 2;
      m_FBar = zeros(ntens,1);
   case 2   %Plane stress ó en el caso que se considere modificar también la Fzz del tensor.
      J0 = Fzz*det(m_F2D);
      pot = 3;
   otherwise
      error('Matrices Elementales FBar_q1: Tipo de cálculo del FBar no implementado.')
end
%Determinación de la inversa de F0
%Inversa de la matriz 2D (F02D)
m_InvF0 = m_F2D\eye(2,2);
%Si no se produce deformación en Fzz viene la identidad.
InvF0zz = 1/Fzz;
%En notación de Voigt, ordenada como [InvF0xx,InvF0yy,InvF0zz,InvF0yx,InvF0xy]. Notar que se invierte los
%términos de corte porque así se necesita en el cálculo de la matriz tangente elemental.
m_InvF0 = [m_InvF0(1,1);m_InvF0(2,2);InvF0zz;m_InvF0(2,1);m_InvF0(1,2)];
%Producto de la matriz de deformación con inv(F0)
m_BInvF0 = B'*m_InvF0;

% Puntos de gauss de 1 a 4. 
for iPG = 1:nPG
   
   e_VG.iPG = iPG;

   B = m_Be(:,:,iPG);
   
   eps_fluct(:,iPG) = B*u;

   % Gradiente de deformación aplicado a cada punto de Gauss.
   m_F = DefMacro+eps_fluct(:,iPG);
   eps_new(:,iPG) = m_F;
   
   %Matriz 2D de F
   m_F2D = [m_F(1),m_F(4);m_F(5),m_F(2)];
   Fzz = m_F(3);
%    tipoFBar=2;
   switch tipoFBar
      case 1  %Plane deformation
         %Se utiliza lo propuesto por Souza, de no modificar la componente Fzz del tensor.         
         % Determinante del gradiente de deformación.
         J = det(m_F2D);
         %Gradiente de deformación FBar
         facFBar = (J0/J)^(1/pot);
         m_FBar([1;2;4;5]) = facFBar*m_F([1;2;4;5]);
         %Por si acaso se mantiene el valor que componente zz, aunque debería ser 1.
         m_FBar(3) = Fzz;
      case 2   %Plane stress ó en el caso que se considere modificar también la Fzz del tensor.
         % Determinante del gradiente de deformación.
         J = Fzz*det(m_F2D);            
         % Gradiente de deformación FBar
         facFBar = (J0/J)^(1/pot);
         m_FBar = facFBar*m_F;
      otherwise
         error('Matrices Elementales FBar_q1: Tipo de cálculo del FBar no implementado.')
   end

   % Modelo constitutivo
   switch conshyp
      case {50,55}  %Modelo multiescala clásico
         [ct,sigma_new(:,iPG),hvar_new(:,iPG)] = f_RMap_ME(m_FBar,hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 52  %Elasticidad usando el tensor elástico obtenido de una homogenización de una microcelda.
         [ct,sigma_new(:,iPG)] = rmap_elast(m_FBar,e_DatMatSet.ce);
      case 100 %Elastic Material neo-Hookean
         [ct,sigma_new(:,iPG)] = f_RMap_NeoHook_DefPl(m_FBar,e_DatMatSet,e_VG);
         f_Const = @(m_F)f_RMap_NeoHook_DefPl(m_F,e_DatMatSet,e_VG);
      case 110 %Large deformations J2 Plasticity
         [ct,sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG),strsg,dmatx] = ...
            f_RMapPlastJ2LD(m_FBar,hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
         %Tensor tangente elástico de Henky en función de G y K.
%          G=e_DatMatSet.young/2/(1+e_DatMatSet.poiss); 
%          K=e_DatMatSet.young/3/(1-2*e_DatMatSet.poiss);
%          d1=4/3*G+K; d2=K-2/3*G; d3=G;
%          c=[d1,d2,0,d2;d2,d1,0,d2;0,0,d3,0;d2,d2,0,d1];
%          f_Const = @(m_FBar,m_F)f_RMapPlastJ2LD(m_FBar,hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      otherwise
         error('Matrices Elementales FBar_q1: Modelo constitutivo no definido.')
   end

   %%
   %Determinación numérica del tensor tangente constitutivo
%    dF = 1e-8;
%    m_DPF = zeros(ntens,ntens);
%    m_Fdf = m_FBar;
%    for i = 1:ntens
%       m_Fdf(i) = m_FBar(i)+dF;
%       [~,m_P1] = f_Const(m_Fdf);
%       m_Fdf(i) = m_FBar(i)-dF;
%       [~,m_P0] = f_Const(m_Fdf);
%       %Diferencia finita centrada
%       m_DPF(:,i) = (m_P1-m_P0)/2/dF;
%       %Se recupera el valor original para realizar los siguientes incrementos.
%       m_Fdf(i) = m_FBar(i);
%    end
%    ct = m_DPF;
   %ct = zeros(ntens,ntens);
   
   %%
%    %Determinación numérica del tensor D[Tau]/D[grads[u]]
%    %Determinación de la inversa de F0
%    %Inversa de la matriz 2D (F02D)
%    m_InvF = m_F2D\eye(2,2);
%    %Si no se produce deformación en Fzz viene la identidad.
%    InvFzz = 1/Fzz;
%    %En notación de Voigt, ordenada como [InvF0xx,InvF0yy,InvF0zz,InvF0yx,InvF0xy]. Notar que se invierte los
%    %términos de corte porque así se necesita en el cálculo de la matriz tangente elemental.
%    m_InvF = [m_InvF(1,1);m_InvF(2,2);InvFzz;m_InvF(2,1);m_InvF(1,2)];
%    m_gradu = reshape(u,e_VG.ndime,[])*reshape(B'*m_InvF,e_VG.ndime,[])';
%    m_IdMenosgradu = [1-m_gradu(1,1);1-m_gradu(2,2);1-0;-m_gradu(1,2);-m_gradu(2,1)];
%    m_IncFIncgradu = m_IdMenosgradu*m_F';
%    %Se divide por dos los términos de corte de m_IncFIncgradu para considerar que se está perturbando un
%    %vector de Voigt, y por lo tanto los términos de corte están multiplicados por 2.
%    m_IncFIncgradu = [m_IncFIncgradu(1,1),m_IncFIncgradu(4,5),(m_IncFIncgradu(1,5)+m_IncFIncgradu(4,1))/2,0;
%       m_IncFIncgradu(5,4),m_IncFIncgradu(2,2),(m_IncFIncgradu(5,2)+m_IncFIncgradu(2,4))/2,0;
%       0,0,0,m_IncFIncgradu(3,3);
%       m_IncFIncgradu(1,4),m_IncFIncgradu(4,2),(m_IncFIncgradu(1,2)+m_IncFIncgradu(4,4))/2,0;
%       m_IncFIncgradu(5,1),m_IncFIncgradu(2,5),(m_IncFIncgradu(5,5)+m_IncFIncgradu(2,1))/2,0];
%    %
%    %
%    dGSu = 1e-8;
%    m_DTGSu = zeros(ntens-1,ntens-1);
%    m_Incgradu = zeros(ntens-1,1);
%    for i = 1:ntens-1
%       m_Incgradu(i) = dGSu;
%       m_IncF = m_IncFIncgradu*m_Incgradu;
%       %m_Fdf = m_FBar+m_IncF;
%       m_Fdf = m_F+m_IncF;
% %       m_F2D = [m_Fdf(1),m_Fdf(4);m_Fdf(5),m_Fdf(2)];
% %       Fzz = m_Fdf(3);
% %       % Determinante del gradiente de deformación.
% %       J = Fzz*det(m_F2D);            
% %       % Gradiente de deformación FBar
% %       facFBar = (J0/J)^(1/pot);
% %       m_Fdf = facFBar*m_F;
%       [~,~,~,~,m_T1] = f_Const(m_Fdf);
%       %m_Fdf = m_FBar-m_IncF;
%       m_Fdf = m_F-m_IncF;
% %       m_F2D = [m_Fdf(1),m_Fdf(4);m_Fdf(5),m_Fdf(2)];
% %       Fzz = m_Fdf(3);
% %       % Determinante del gradiente de deformación.
% %       J = Fzz*det(m_F2D);            
% %       % Gradiente de deformación FBar
% %       facFBar = (J0/J)^(1/pot);
% %       m_Fdf = facFBar*m_F;
%       [~,~,~,~,m_T0] = f_Const(m_Fdf);
%       %Diferencia finita centrada
%       m_DTGSu(:,i) = (m_T1-m_T0)/2/dGSu;
%       %Se recupera el valor original para realizar los siguientes incrementos.
%       m_Incgradu(i) = 0;
%    end
%    %Tensor D[P]/D[F], pero usando el c calculado a partir de perturbarción.
%    %Inversa de la matriz 2D (F02D)
%    m_F2DBar = [m_FBar(1),m_FBar(4);m_FBar(5),m_FBar(2)];
%    FzzBar = m_FBar(3);
%    m_InvFBar = m_F2DBar\eye(2,2);
%    %m_InvF = m_F2D\eye(2,2);
%    %Si no se produce deformación en Fzz viene la identidad.
%    InvFzzBar = 1/FzzBar;
%    %InvFzz = 1/Fzz;
%    %
%    m_DPF2 = f_DTau_DgraduaDP_DF(strsg,m_DTGSu,m_InvFBar,InvFzzBar);
%    %ct = m_DPF2;
   
   %%
   %Se almacena para los tensor tangente constitutivos para realizar homogeneización y análisis de
   %bifurcación
   if esImplex
      % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
      %tangentes implícitas para el análisis de bifurcación.
      %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando está seleccionado
      %el implex.
      m_TensorTang(:,:,iPG) = ct.Implex;
      %Se almacena los tensores tangentes constitutivo implícitos para análisis de bifurcación como si fuera
      %PG adicionales, tantos como nPG. Se almacena en los índices (:,:,nPG+1:2*nPG).
      %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
      %tercera dimensión de m_TensorTang (size(m_TensorTang,3)).
      m_TensorTang(:,:,iPG+nPG) = ct.Impli;
      %En los cálculos para el ensamblaje se utiliza el implex.
      ct = ct.Implex;
   else
      if isfield(ct,'Impli')  
         ct = ct.Impli;
      end
      m_TensorTang(:,:,iPG) = ct;
   end
   
   % Cálculo de fint
   %delta_sigma = sigma_new(:,iPG)-sigma_old(:,iPG);   
   %fint = fint+B'*delta_sigma*m_pesoPG(iPG);
   fint = fint+B'*sigma_new(:,iPG)*m_pesoPG(iPG);

   %Determinación de la inversa de F0
   %Inversa de la matriz 2D (F02D)
   m_InvF = m_F2D\eye(2,2);
   %Si no se produce deformación en Fzz viene la identidad.
   InvFzz = 1/Fzz;
   %En notación de Voigt, ordenada como [InvF0xx,InvF0yy,InvF0zz,InvF0yx,InvF0xy]. Notar que se invierte los
   %términos de corte porque así se necesita en el cálculo de la matriz tangente elemental.
   m_InvF = [m_InvF(1,1);m_InvF(2,2);InvFzz;m_InvF(2,1);m_InvF(1,2)];
   %Producto de la matriz de deformación con inv(F0)
   m_BInvF = B'*m_InvF;
   
   % Cálculo de Ke
   switch tipoFBar
      case 1  %Plane deformation
         %Se utiliza lo propuesto por Souza, de no modificar la componente Fzz del tensor.  
         %Acá se está asumiendo, como es normal para deformación plana, que B, en zz tiene ceros.
         kt = kt+facFBar*B'*(ct(:,[1;2;4;5])*m_F([1;2;4;5])*(m_BInvF0-m_BInvF)'/pot+...
            ct*B)*m_pesoPG(iPG);
      case 2   %Plane stress ó en el caso que se considere modificar también la Fzz del tensor.
         kt = kt+facFBar*B'*(ct*m_F*(m_BInvF0-m_BInvF)'/pot+ct*B)*m_pesoPG(iPG);
         %kt = kt+facFBar*B'*ct*(m_F*(m_BInvF0-m_BInvF)'/pot+B)*m_pesoPG(iPG);
      otherwise
         error('Matrices Elementales FBar_q1: Tipo de cálculo del FBar no implementado.')
   end
   
   
   
   
   
   % Calculo de la energia interna
   switch conshyp   %JLM habilite esto que estaba comentado
%        case 1 % Elastic material
%            IntEnergyNew(iPG) = 0;
       case 2 % Plastic material W = 0.5*(e_e*C*e_e)+0.5*K*chi^2
           IntEnergyNew(iPG) = 0.5*((eps_new(:,iPG)-hvar_new(1:4,iPG))'*e_DatMatSet.ce*(eps_new(:,iPG)-hvar_new(1:4,iPG))) ...
               +0.5*hvar_new(12,iPG)*(hvar_new(5,iPG))^2;
       case 110 
%%            
%            Free Energy Psi = Psi^e + Psi^p
%                           Psi^e = 1/2*eps:D:eps ; con eps = 1/2*ln(B), B=F*F', D=tensor elastico de Hencky
%                           Psi^p = Sy*alpha + 1/2*Thita*H*alpha^2 + (K_oo
%                           - K_o)*(alpha + 1/delta*exp(-delta*alpha))
%                        dPsi^p / d alpha = K(alpha) de Simo !!!

%          Tensor tangente elástico de Henky en función de G y K.
%          G = e_DatMatSet.young/2/(1+e_DatMatSet.poiss); 
%          K = e_DatMatSet.young/3/(1-2*e_DatMatSet.poiss);
%          d1= 4/3*G+K; d2=K-2/3*G; d3=G;
%          D = [d1 d2 0 d2 d2; d2 d1 0 d2 d2; 0 0 d3 0 0;  d2 d2 0 d1 d2; d2 d2 0 d2 d1];
%          
%          eps = eps_new(:,iPG) ;
%          alpha = hvar_new(1,iPG) ;
%          F = [eps(1) eps(4); eps(5) eps(2) ] ;
%          B = F*F'  ;  
%          [v,a]=eig(B) ;
%          ln_B = log(a(1,1))*v(:,1)*v(:,1)'+log(a(2,2))*v(:,2)*v(:,2)';
% %          B = [eps(1)^2+eps(4)^2; eps(2)^2+eps(5)^2; eps(3)^2; ...
% %              eps(1)*eps(5)+eps(2)*eps(4); eps(1)*eps(5)+eps(2)*eps(4) ];
% %          eps_ln = real(0.5*log(B)) ;% MAL el log no se hace asi
% %         B=F*F'  ;  [v,a]=eig(B) ; B=a(1,1)*v(:,1)*v(:,1)'+a(2,2)*v(:,2)*v(:,2)';
%          eps_ln = 0.5*[ln_B(1,1) ln_B(2,2) 0 ln_B(1,2) ln_B(2,1)]';
%          Psi_e = 0.5*eps_ln'*D*eps_ln ;
           
%        Elastic energies
        Psi_e     = hvar_new(7,iPG) ;
        Psi_e_vol = hvar_new(8,iPG) ;
        Psi_e_dev = hvar_new(9,iPG) ;
             
%La función de módulo de endurecimiento isotrópico es (pág. 91 del Computational Elasticity - Simo):
         %K(alpha) = YIELD+TIT*HBA*alfa+(KIN-KCE)*(1-exp(-DEL*alfa))
         alpha = hvar_new(1,iPG) ;
         YIELD = e_DatMatSet.yield ;
         TIT   = e_DatMatSet.tit ;
         HBA   = e_DatMatSet.hba ;  
         KIN   = e_DatMatSet.kin ;
         KCE   = e_DatMatSet.kce ;
         DEL   = e_DatMatSet.del ;
         
         Psi_p = YIELD*alpha + 0.5*TIT*HBA*alpha^2 + (KIN-KCE)*(alpha+1/DEL*exp(-DEL*alpha)) ;
         if alpha ~= 0
             Psi_p=Psi_p;
         end
         
%        Total Energy         
         IntEnergyNew((iPG-1)*4+1) = Psi_e ;
         IntEnergyNew((iPG-1)*4+2) = Psi_e_vol ; 
         IntEnergyNew((iPG-1)*4+3) = Psi_e_dev ;
         IntEnergyNew((iPG-1)*4+4) = Psi_p ;
          
         
%            IntEnergyNew(iPG) = 0.5*(sigma_new(:,iPG)'*eps_new(:,iPG));  % energia interna
%            IntEnergyNew(iPG) = 0.5*(sigma_new(:,iPG)'*[eps_new(:,iPG)-[1 1 1 0 0]']);  % energia interna
%            IntEnergyNew(iPG) = 0.5*([eps_new(:,iPG)-[1 1 1 0 0]']'*ct*[eps_new(:,iPG)-[1 1 1 0 0]']);  % energia interna           
%            IntEnergyNew(iPG) = 0.5*((eps_new(:,iPG)-hvar_new(1:5,iPG))'*e_DatMatSet.ce*(eps_new(:,iPG)-hvar_new(1:4,iPG))) ...
%                +0.5*hvar_new(12,iPG)*(hvar_new(5,iPG))^2;
%        case 11  % Continuum damage model
%            IntEnergyNew(iPG) = 0.5*((1-aux_var(1,iPG))*eps_new(:,iPG)'*...
%                e_DatMatSet.ce*eps_new(:,iPG));  % energia interna
%        case 50 % Multiscale classic model
%            
%            %IntEnergyNew(iPG) = 0.5*(sigma_new(:,iPG)'*eps_new(:,iPG));  % energia interna
%            
%            % Flag to split the snapshots (elastic & inelastic part)
%            IntEnergyNew(iPG) = max(max(hvar_new(1).e_VarEst(1).hvar(13:e_DatMatSet.e_DatSet(1).e_DatMat.sihvarpg:end,:)));
%        otherwise
%            IntEnergyNew(iPG) = 0.5*((1-aux_var(1,iPG))*eps_new(:,iPG)'*...
%                e_DatMatSet.ce*eps_new(:,iPG));  % energia interna
   end
   
   
   
end

% Se ordena las matrices como vectores columnas
sigma_new = sigma_new(:);
hvar_new = hvar_new(:);
aux_var = aux_var(:);
eps_new = eps_new(:);
eps_fluct = eps_fluct(:);
   
end