function [m_kt,fint,sigma_new,eps_new,eps_fluct,hvar_new,m_VarHistElemNew,m_VarAuxPG,m_VarAuxElem,m_TensorTang,...
   m_resb,m_kbu,m_invkbb] = ...
   f_MatElem_SDA_tria_t1(...
   u,m_Beta,hvar_old,m_VarHistElemOld,m_VarAuxPG,m_VarAuxElem,e_DatElemSet,e_DatMatSet,m_Be,m_DetJe,DefMacro,sigma_old,...
   eps_old,e_VG)
   
   %Solo se considera el caso de dos puntos de gauss, uno para la parte continua (regular) y otro
   %para la parte discontinua (singular).
   %Ver si despu�s no hacer que sea posible ingresar la cantidad de grados de libertad que se quiera
   %antes de la bifurcaci�n y se quede con dos en el centro despu�s de la bifurcaci�n (�sirve para
   %algo en el tri�ngulo lineal?).
   
   %condBif = 0: No se produjo todav�a la bifurcaci�n.
   %condBif = 1: Se produjo la bifurcaci�n pero todav�a no se activa la cinem�tica de
   %discontinuidades fuertes y el modelo cohesivo.
   %condBif = 2: Se activa la discontinuidad fuerte y modelo cohesivo despu�s de la bifurcaci�n. A
   %partir de este punto se empieza a medir el salto, con valor inicial nulo.
   %condBif = 3: Cuando ya se pasa un paso desde la activaci�n de la discontinuidad fuerte y el modelo
   %cohesivo. Sirve para realizar operaciones una sola vez despu�s que se activ� la bifurcaci�n.
   
   % Recuperaci�n de variables
   ntens = e_VG.ntens;
   ndime = e_VG.ndime;
   %
   dofpe = e_DatElemSet.dofpe;
   npg = e_DatElemSet.npg; 
   wg = e_DatElemSet.wg;
   anchoLoc = e_DatElemSet.anchoLoc;
   %
   sihvarpg = e_DatMatSet.sihvarpg;
   siavarpg = e_DatMatSet.siavarpg;
   conshyp  = e_DatMatSet.conshyp;
   aux_var = reshape(m_VarAuxPG,siavarpg,npg);
   esImplex = e_DatMatSet.esImplex;
   %
   condBif = m_VarAuxElem(1);

   % Inicializaciones
   sigma_new = zeros(ntens,npg);
   eps_new = zeros(ntens,npg);
   eps_fluct = zeros(ntens,npg);
   %hvar_new = zeros(sihvarpg,npg);
   hvar_new = f_DefinicionhVar(conshyp,sihvarpg,npg);
   if esImplex
      m_TensorTang = zeros(ntens,ntens,2*npg);
   else
      m_TensorTang = zeros(ntens,ntens,npg);
   end
   
   % Redimensionado de matrices
   hvar_old = reshape(hvar_old,sihvarpg,npg);
   %hvar_new = reshape(hvar_new,sihvarpg,npg);
   sigma_old = reshape(sigma_old,ntens,npg);
   eps_old = reshape(eps_old,ntens,npg);
   m_VarAuxPG = reshape(m_VarAuxPG,siavarpg,npg);
 
   %Recuperaci�n de variables
   m_Normal = reshape(m_VarAuxElem(2:9),ntens,ndime);  

   %En el peso de gauss ya viene multiplicado el espesor.
   m_pesoPG = m_DetJe.*wg;
   %Como se tiene dos puntos de gauss en la misma posici�n, para calcular el volumen hay que considerar que
   %los dos puntos de gauss est�n en la misma posici�n, y por lo tanto hay que usar uno solo.
   %volElem = sum(m_pesoPG);
   %volElem = m_pesoPG(1);
   
   %Se asume que la matriz de deformaci�n precalculada es la misma para ambos PG, al estar en la
   %misma posici�n.
   B = m_Be(:,:,1);
   %Esta parte de la deformaci�n es la misma para ambos PG ya que la matriz de deformaci�n es la
   %misma.
   m_EpsNewR = B*u;
   %Habr�a que ver en este caso c�mo queda definida la deformaci�n regular (esto servir�a solo para 
   %el caso que usa SDA a nivel micro).
   %eps_fluct(:,1) = m_EpsNewR;
   %Normalmente, si no se llama un problema multiescala concatenado uno dentro del otro, a este
   %nivel, DefMacro es a nivel macro, por lo tanto nula.
   m_EpsNewR = DefMacro(:,1)+m_EpsNewR;
   %eps_new(:,1) = m_EpsNewR;
   if condBif<2
      noActBif = m_VarAuxElem(21);
      if condBif==1&&noActBif==1
         %A los elementos con noActBif==1, despu�s que bifurc� se impone que sean el�stico, as� forzar su
         %descarga. Son elementos en la zona de exclusi�n.
         e_VG.elast = 1;
      end         
      switch conshyp
         case 10   %Da�o isotr�pico regularizado
            [ctR,sigma_new(:,1),eps_new(:,1),hvar_new(:,1),m_VarAuxPG(:,1)] = rmap_damage_reg(...
               m_EpsNewR,hvar_old(:,1),e_DatMatSet,e_VG);
            dsigmaR = sigma_new(:,1)-sigma_old(:,1);
         case 12   %Da�o isotr�pico s�lo tracci�n regularizado            
            [ctR,sigma_new(:,1),eps_new(:,1),hvar_new(:,1),m_VarAuxPG(:,1)] = ...
               RMapDanoSTraccReg(...
               m_EpsNewR,hvar_old(:,1),aux_var(:,iPG),e_DatMatSet,e_VG);
            dsigmaR = sigma_new(:,1)-sigma_old(:,1);
         case 51   %Modelo multiescala cohesivo
            %fprintf('**-- Inicio del return mapping del modelo multiescala\n')
%             cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
%             listeners = cmdWinDoc.getDocumentListeners;
%             jFxCommandArea = listeners(3);
%             set(jFxCommandArea,'Background','red');
            %
            e_VG.iPG = 1;
            m_IDefR = m_EpsNewR-eps_old(:,1);
            m_IDefS = zeros(ntens,1);
            [ctR,sigma_new(:,1),hvar_new(:,1),m_VarAuxPG(:,1)] = ...
              f_RMap_MECohesivo(...
               m_IDefR,m_IDefS,sigma_old(:,1),hvar_old(:,1),m_VarAuxPG(:,1),e_DatMatSet,condBif,e_VG);
            eps_new(:,1) = m_EpsNewR;
            %
%             set(jFxCommandArea,'Background','yellow');
            %fprintf('**-- Fin del return mapping del modelo multiescala\n')
          otherwise
            error(['Matrices Elementales SDA Tria_t1: Condici�n estable del material: ',...
               'Modelo constitutivo no definido.'])
      end
      %Se almacena para los tensor tangente constitutivos para realizar homogeneizaci�n y an�lisis de
      %bifurcaci�n
      if esImplex
         % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
         %tangentes impl�citas para el an�lisis de bifurcaci�n.
         %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando est�
         %seleccionado el implex, a�n cuando despu�s de la bifurcaci�n ya no necesita los tensores
         %constitutivos impl�citos.
         %Se guarda el tensor constitutivo singular como el mismo que el regular (ctS = ctR) 
         m_TensorTang(:,:,1:2) = repmat(ctR.Implex,[1,1,2]);
         %Se almacena los tensores tangentes constitutivo impl�citos para an�lisis de bifurcaci�n como si
         %fuera PG adicionales, tantos como nPG. Se almacena en los �ndices (:,:,nPG+1:2*nPG).
         %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
         %tercera dimensi�n de m_TensorTang (size(m_TensorTang,3)).
         m_TensorTang(:,:,3:4) = repmat(ctR.Impli,[1,1,2]);
         %En los c�lculos para el ensamblaje se utiliza el implex.
         ctR = ctR.Implex;
      else
         %Se guarda el tensor constitutivo singular como el mismo que el regular (ctS = ctR)
         m_TensorTang(:,:,1) = ctR;
         m_TensorTang(:,:,2) = ctR;
      end
      %
    %  fint = B'*dsigmaR*m_pesoPG(1);
      fint = B'*sigma_new(:,1)*m_pesoPG(1);
      m_kt = B'*ctR*B*m_pesoPG(1);
      m_resb = zeros(ndime,1);
      m_kbu = zeros(ndime,dofpe);
      m_invkbb = zeros(ndime,ndime);
     %  m_TraccNew = zeros(ndime,1);   ### alfredo
      m_TraccNew = m_Normal'*sigma_new(:,1);
   elseif condBif>1
      %Recuperaci�n de variables
     %  m_Normal = reshape(m_VarAuxElem(2:9),ntens,ndime);  
      m_GradPhi = reshape(m_VarAuxElem(10:17),ntens,ndime);
      %Deformaci�n regular
      m_EpsNewR = m_EpsNewR-m_GradPhi*m_Beta;
      switch conshyp
         case 10    %Da�o isotr�pico regularizado
            % Punto de gauss regular
            %Para obtener una longitud del elemento se considera como un prisma recto de altura 1
            %con bases de tri�ngulos rect�ngulos.
            %hElemReg = sqrt(2*volElem/1);
            %[ctR,sigma_new(:,1),eps_new(:,1),hvar_new(:,1),m_VarAuxPG(:,1)] = rmap_damage_reg(...
            %   m_EpsNewR,hvar_old(:,1),e_DatMatSet,hElemReg,e_VG);
            %Se fuerza descarga el�stica
            ctR = e_DatMatSet.ce*hvar_old(2,1)/hvar_old(1,1);
            sigma_new(:,1) = ctR*m_EpsNewR;
            hvar_new(:,1) = hvar_old(:,1);
            eps_new(:,1) = m_EpsNewR;
            %fload (descarga el�stica)
            hvar_new(4,1) = -1;
            dsigmaR = sigma_new(:,1)-sigma_old(:,1);
            % Punto de gauss singular
            m_EpsNewS = m_EpsNewR+m_Normal*m_Beta/anchoLoc;  
            %Par�metro de ablandamiento discreto (se pasa con la longitud de regularizaci�n)
            e_DatMatSet.m_hReg(e_VG.iElemSet) = anchoLoc; 
            [ctS,sigma_new(:,2),eps_new(:,2),hvar_new(:,2),m_VarAuxPG(:,2)] = rmap_damage_reg(...
               m_EpsNewS,hvar_old(:,2),e_DatMatSet,e_VG);
            dsigmaS = sigma_new(:,2)-sigma_old(:,2);
         case 12    %Da�o isotr�pico s�lo tracci�n regularizado
            % Punto de gauss regular
            %e_VG.iPG = 1;
            %Se exige que los modelos constitutivos micros se comporten el�sticamente.
            e_VG.elast = 1;
            [ctR,sigma_new(:,1),eps_new(:,1),hvar_new(:,1),m_VarAuxPG(:,1)] = ...
               RMapDanoSTraccReg(...
               m_EpsNewR,hvar_old(:,1),aux_var(:,iPG),e_DatMatSet,e_VG);
            %Se fuerza descarga el�stica
            dsigmaR = sigma_new(:,1)-sigma_old(:,1);
            % Punto de gauss singular
            %e_VG.iPG = 2;
            e_VG.elast = 0;
            m_EpsNewS = m_EpsNewR+m_Normal*m_Beta/anchoLoc;  
            %Par�metro de ablandamiento discreto (se pasa con la longitud de regularizaci�n)
            e_DatMatSet.m_hReg(e_VG.iElemSet) = anchoLoc; 
            [ctS,sigma_new(:,2),eps_new(:,2),hvar_new(:,2),m_VarAuxPG(:,2)] = ...
               RMapDanoSTraccReg(...
               m_EpsNewS,hvar_old(:,2),aux_var(:,iPG),e_DatMatSet,e_VG);
            dsigmaS = sigma_new(:,2)-sigma_old(:,2);
         case 51   %Modelo multiescala cohesivo
            %fprintf('**-- Inicio del return mapping del modelo multiescala\n')
%             cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
%             listeners = cmdWinDoc.getDocumentListeners;
%             jFxCommandArea = listeners(3);
%             set(jFxCommandArea,'Background','red');
            %
            % Punto de gauss regular ######################
            e_VG.iPG = 1;
            e_VG.elast = 1;  % respuesta elastica
            condBifR = 0;
            %Incremento de la deformaci�n regular
            m_IDefR = m_EpsNewR-eps_old(:,1);
            m_IDefS = zeros(ntens,1);
            % Para cambiar la energ�a de fractura del PG regular
            [ctR,sigma_new(:,1),hvar_new(:,1),m_VarAuxPG(:,1)] = f_RMap_MECohesivo(...
               m_IDefR,m_IDefS,sigma_old(:,1),hvar_old(:,1),m_VarAuxPG(:,1),e_DatMatSet,condBifR,e_VG);
            %
            eps_new(:,1) = m_EpsNewR;
            % Punto de gauss singular ######################
            e_VG.iPG = 2;
            e_VG.elast = 0;    % respuesta danio
            condBifS = condBif;
            %Por ahora hasta que se programe bien el tensor tangente, se utiliza como ancho localizado la
            %media de las longitudes micro.
            anchoLoc = hvar_old(:,2).lMicro;
            %Incremento de la parte singular (b*n/k) de la deformaci�n
            %m_IDefS = m_Normal*m_Beta/anchoLoc-(eps_old(:,2)-eps_old(:,1));
            m_BetaOld = m_VarHistElemOld(1:ndime);
            m_IDefS = m_Beta-m_BetaOld;
            %Funci�n del modelo constitutivo multiescala cohesivo
            [ctS,m_TraccNew,hvar_new(:,2),m_VarAuxPG(:,2)] = ...
             f_RMap_MECohesivo(...
                m_IDefR,m_IDefS,sigma_old(:,2),hvar_old(:,2),m_VarAuxPG(:,2),e_DatMatSet,condBifS,e_VG);
%             set(jFxCommandArea,'Background','yellow');
         otherwise
            error('Matrices Elementales SDA Tria_t1: Modelo constitutivo no definido.')
      end
      %Se almacena para los tensor tangente constitutivos para realizar homogeneizaci�n y an�lisis de
      %bifurcaci�n
      if esImplex
         % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
         %tangentes impl�citas para el an�lisis de bifurcaci�n.
         %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando est�
         %seleccionado el implex.
         m_TensorTang(:,:,1) = ctR.Implex;
         m_TensorTang(:,:,2) = ctS.Implex;
         %Se almacena los tensores tangentes constitutivo impl�citos para an�lisis de bifurcaci�n como si
         %fuera PG adicionales, tantos como nPG. Se almacena en los �ndices (:,:,nPG+1:2*nPG).
         %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
         %tercera dimensi�n de m_TensorTang (size(m_TensorTang,3)).
         m_TensorTang(:,:,3) = ctR.Impli;
         m_TensorTang(:,:,4) = ctS.Impli;
         %En los c�lculos para el ensamblaje se utiliza el implex.
         ctR = ctR.Implex;
         ctS = ctS.Implex;
      else
         % Se almacena para la homogenizaci�n del tensor tangente
         m_TensorTang(:,:,1) = ctR;
         %m_TensorTang(:,:,2) = ctS;
      end
      %Se asume que la normal no var�a en el tiempo desde que se impone el SDA y por lo tanto las
      %tracciones siempre son T = sigma*n, y no T = T_old+dsigma*n (sino hay que llevar como
      %variable hist�rica la tracci�n).
      %Se usa las tensiones del punto de gauss singular, que hasta que converja el problema es
      %distinta al regular, pero si se utiliza posterior a que esto suceda las traciones obtenidas
      %con ambos PG es la misma.
      %m_VarAuxElem(19:20) = m_Normal'*sigma_new(:,2);
      %Se recupera la variable de la tracci�n, ya que ahora el modelo constitutivo devuelve directamente la
      %tracci�n.
  %    m_TraccOld = m_VarHistElemOld(ndime+1:2*ndime);
      %Se la sigue almacenando en las variables auxiliares para que lo imprima correctamente.
      %m_VarAuxElem(19:20) = m_TraccOld+dsigmaS;
  %    m_TraccNew = m_TraccOld+dsigmaS;
      
     % fu = B'*dsigmaR*m_pesoPG(1); ###alfredo
      fu = B'*sigma_new(:,1)*m_pesoPG(1);
      
     % fb = (dsigmaS-m_Normal'*dsigmaR)*m_pesoPG(1);
      fb = (m_TraccNew-m_Normal'*sigma_new(:,1))*m_pesoPG(1);
      m_kuu = B'*ctR*B*m_pesoPG(1);
      m_kub = -B'*ctR*m_GradPhi*m_pesoPG(1);
      
%%%%%%%%%%%%%%%%%%%%%%%%%      
%alf      m_Dctn = m_Normal'*(ctS-ctR);
    m_Dctn = ctS(ndime+1:2*ndime,:) +   m_Normal'*( -ctR);
 %%%%%%%%%%%%%%%%%%%%%%%%%      
     
      
      %m_Dctn = (ctS-m_Normal'*ctR);
      %m_Dctn = -m_Normal'*ctR;
      %m_kbu = m_Dctn*B;
      m_kbu = m_Dctn *B*m_pesoPG(1);
      %m_invkbb = inv(m_Normal'*ctS*m_Normal-m_Dctn*m_GradPhi);
      %m_invkbb = (m_Normal'*ctS*m_Normal-m_Dctn*m_GradPhi)\eye(ndime);
      
%%%%%%%%%%%%%%%%%%%%%%%%%      
%alf      m_invkbb = ((m_Normal'*ctS/anchoLoc*m_Normal-m_Dctn*m_GradPhi)*m_pesoPG(1))\eye(ndime);
      m_invkbb = ( (ctS(1:ndime, 1:ndime)- m_Dctn*m_GradPhi)*m_pesoPG(1))\eye(ndime);
%%%%%%%%%%%%%%%%%%%%%%%%%      
      
      
      %m_invkbb = ((ctS/anchoLoc*m_Normal-m_Dctn*m_GradPhi)*m_pesoPG(1))\eye(ndime);
      %El residuo de beta es fb-fextb, pero se asume que fextb directamente es nulo.
      m_resb = fb;
      %
      m_k1 = m_kub*m_invkbb;
      %Fuerza interna y matriz de rigidez elemental del problema condensado (se condesa los grados
      %de libertad internos de salto)
      fint = fu-m_k1*fb;
      m_kt = m_kuu-m_k1*m_kbu;
      %
   end
   
   % Se ordenan las matrices como vectores columnas
   sigma_new = sigma_new(:);
   hvar_new = hvar_new(:);
   m_VarAuxPG = m_VarAuxPG(:);
   eps_new = eps_new(:);
   eps_fluct = eps_fluct(:);
   m_VarHistElemNew(1:ndime) = m_Beta;
   m_VarHistElemNew(ndime+1:2*ndime) = m_TraccNew;
   
end