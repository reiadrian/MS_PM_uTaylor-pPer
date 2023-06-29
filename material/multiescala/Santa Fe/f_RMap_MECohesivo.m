function [m_CTHomog,m_SigmaHomog,hvar_newMacro,aux_varMacro] = ...
  f_RMap_MECohesivo(...
   m_IDefMacroReg,m_IDefMacroLoc,m_sigmaHomog_old,hvar_oldMacro,aux_varMacro,...
   e_DatMatSetMacro,condBif,e_VGMacro)

   %Se recupera variables micro
   xx = e_DatMatSetMacro.xx;
   omegaMicro = e_DatMatSetMacro.omegaMicro;
   e_DatSet = e_DatMatSetMacro.e_DatSet;
   e_VG = e_DatMatSetMacro.e_VG;
   esImplexMacro = e_DatMatSetMacro.esImplex;
   e_VarEst_old = hvar_oldMacro.e_VarEst;
   u = hvar_oldMacro.u;
   c_GdlCond = hvar_oldMacro.c_GdlCond;
   Fint = hvar_oldMacro.Fint;
   m_LinCond = hvar_oldMacro.m_LinCond;
   doff = hvar_oldMacro.doff;
   dofl = hvar_oldMacro.dofl;
   %vfix = hvar_oldMacro.vfix;
   e_VarAux = hvar_oldMacro.e_VarAux;
   c_DefMacro = hvar_oldMacro.c_DefMacro;
   omegaMicroL = hvar_oldMacro.omegaMicroL;
   m_ElemLoc = hvar_oldMacro.m_ElemLoc;
   lMacro = hvar_oldMacro.lMacro;
   lMicro = hvar_oldMacro.lMicro;
%    m_TensProy = hvar_oldMacro.m_TensProy;
   c_NormalesMicro = hvar_oldMacro.c_NormalesMicro;
   longFis = hvar_oldMacro.longFis;
   facNormMicro = hvar_oldMacro.facNormMicro;


   % VARIABLES GLOBALES
   nElem = e_VG.nElem;
   ntens = e_VG.ntens;
   ndoft = e_VG.ndoft;
   nSet = e_VG.nSet;
   nDime = e_VG.ndime;

   % INICIALIZACION DE VARIABLES
   % Vector de fzas internas
   %Fint = zeros(ndoft,1);   
   % Vector de incrementos de desplazamientos en el paso de tiempo previo
   Du_step_old = zeros(ndoft,1);  
   Fext = zeros(ndoft,1);

   %Por si es necesario el paso tiempo a nivel micro.
   e_VG.istep = e_VGMacro.istep;
   %Se guarda que se quiere que el modelo constitutivo se comporte el�sticamente.
   e_VG.elast = e_VGMacro.elast;
   %Para impresi�n y debug se guarda la iteraci�n macro.
   e_VG.iterMacro = e_VGMacro.iter;
   %Para imprensi�n se guarda el n�mero de elemento macro y n�mero de PG macro.
   e_VG.iElemNumMacro = e_VGMacro.iElemNum;
   e_VG.iPGMacro = e_VGMacro.iPG;
   %Nombre interno de la celda unitaria 
   %(hacer esto en cada iteraci�n puede ser medio lento, ver que hacer sino).
   %e_VG.fileCompleto = [e_VG.fileCompleto,'_EM',int2str(e_VGMacro.iElemNum),'PGM',int2str(e_VGMacro.iPG)];
   %Se guarda los datos del matlabpool (se utiliza para imprimir de pasos y iteraciones a nivel micro que no
   %convergieron).
   e_VG.nLab = e_VGMacro.nLab;
   e_VG.tipoMPool = e_VGMacro.tipoMPool;
   
   % DEFORMACION MACRO TOTAL A APLICAR EN CADA PG  c_DefMacro
   if condBif>1      
      c_DefMacro = cellfun(@(x)bsxfun(@plus,x,m_IDefMacroReg),c_DefMacro,'UniformOutput',false);      
      c_NormalesMicro = hvar_oldMacro.c_NormalesMicro;
      %En m_IDefMacroLoc viene el incremento del salto beta.
      for iSet = 1:e_VG.nSet
         nElem = e_DatSet(iSet).nElem;
         nPG = e_DatSet(iSet).e_DatElem.npg;
         m_NormalesMicro = c_NormalesMicro{iSet};
         m_lMicro = c_NormalesMicro{iSet,2};
         m_DefMacro = c_DefMacro{iSet};
         for iElem = 1:nElem
            for iPG = 1:nPG
               %Se asume que las normales micro que no pertenecen al dominio localizado son nulas.
                m_DefMacro(:,iPG,iElem) = m_DefMacro(:,iPG,iElem)+...
                  m_NormalesMicro(:,:,iElem)*m_IDefMacroLoc/m_lMicro(iElem)/facNormMicro;
            end
         end 
         c_DefMacro{iSet} = m_DefMacro;
      end      
   else
      %Se asume que IDefMacroLocPG es nula previo a la bifurcaci�n y a la imposici�n de la SD.
      c_DefMacro = cellfun(@(x)bsxfun(@plus,x,m_IDefMacroReg),c_DefMacro,'UniformOutput',false);   
   end
   
   % ESQUEMA DE NEWTON-RAPHSON
   %ticIDNewMicr = tic;
   %fprintf('Tiempo del Newton micro: %f\n',toc(ticIDNewMicr));
   Du_step_new = zeros(ndoft,1);
   [u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] =  newton_raphson_MICRO(...
      xx,m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
      c_DefMacro,e_VG);
   %fprintf('Tiempo del Newton micro: %f\n',toc(ticIDNewMicr));
   %Para considerar que en el caso del cuadr�ngulo Q1_Mixed_StrainInjection se debe integrar las tensiones estabilizadas.
   c_Tens_old = cell(nSet,1);
   c_Tens_new = cell(nSet,1);
   for iSet = 1:nSet
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      eltype = e_DatElemSet.eltype;
      switch eltype
         case {2,4,8}
            c_Tens_old{iSet} = e_VarEst_old(iSet).sigma;
            c_Tens_new{iSet} = e_VarEst_new(iSet).sigma;
         case 20
            c_Tens_old{iSet} = e_VarEst_old(iSet).VarHistElem;
            c_Tens_new{iSet} = e_VarEst_new(iSet).VarHistElem;
         otherwise
            error(['Modelo Multiescala: Homogeneizaci�n: Tensiones de homogeneizaci�n: Modelo ',...
               'constitutivo no definido.'])
      end
   end
   
   % INCREMENTO DE TENSIONES HOMEGENEIZADAS
   if condBif>1
      %Esta matriz de los determinante de J con valores no nulos dentro del dominio omegaL podr�a 
      %precalcularse si fuese necesario.
      c_DetJTLoc = arrayfun(@(x)bsxfun(@times,x.m_DetJT,m_ElemLoc(x.m_IndElemSet)),e_DatSet,...
         'UniformOutput',false);
      %Es necesario realizar esta resta en todos los puntos de gauss de los elementos micro porque justo en el
      %paso en que se empieza a homogeneizar en omegaL, la tensi�n previa homogenizada m_sigmaHomog_old, que
      %fue obtenida en una homogeneizaci�n en todo omegaMicro, puede ser distinta a la homogenizaci�n de las
      %tensiones sigma_old en omegaL (al ser dominios distintos). Se debe interpretar que todas las tensiones
      %previas de cada PG de la estructura micro son las tensiones previas del problema macro y no sus valores
      %homogeneizados.
      %c_Dsigma = cellfun(@(x,y)minus(x,y),c_Tens_new,c_Tens_old,'UniformOutput',false);
     
     %%%%%%%%%%%%%%%%%%%
%       c_Tracc=cell(nSet,1);

      m_VarFluc = f_VarFluct(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,m_ElemLoc,e_VG);
      m_CTHomog      = zeros(ntens,ntens);
      m_CTHomog_bb   = zeros(nDime,nDime);
      m_CTHomog_ub   = zeros(nDime,ntens);
      m_IndTens = eye(ntens);
      c_Tracc   = zeros(nDime,1);
      %%%%%%%%%%%%%%%%%%%

      %
      for iSet = 1:e_VG.nSet
         nElem = e_DatSet(iSet).nElem;
         nPG = e_DatSet(iSet).e_DatElem.npg;
         m_NormalesMicro = c_NormalesMicro{iSet,1};
         m_lMicro = c_NormalesMicro{iSet,2};
         m_Tracc = zeros(nDime,nPG,nElem);
         m_Sigma = reshape(c_Tens_new{iSet},ntens,nPG,nElem);
         
         %%%%%%%%%%%%%%%%%%%
         m_CT = c_CT{iSet};
         m_DofElem = e_DatSet(iSet).m_DofElem;
         m_DetJT = e_DatSet(iSet).m_DetJT;
         m_BT = e_DatSet(iSet).m_BT;
         m_CT = c_CT{iSet};
         e_DatElemSet = e_DatSet(iSet).e_DatElem;
         npg = e_DatElemSet.npg;
         wg = e_DatElemSet.wg;
         m_ElemLocSet = m_ElemLoc(e_DatSet(iSet).m_IndElemSet);
         %%%%%%%%%%%%%%%%%
         
         for iElem = 1:nElem
             
             if m_ElemLocSet(iElem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                dofElem = m_DofElem(:,iElem);
                m_pesoPG = m_DetJT(:,iElem).*wg;
                m_VarFlucElem = m_VarFluc(dofElem,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                for iPG = 1:nPG
                     %Se asume que las normales micro que no pertenecen al dominio localizado son nulas y que las
                     %lMicro son Inf.
                     
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     m_Tracc(:,iPG,iElem) = m_NormalesMicro(:,:,iElem)'*m_Sigma(:,iPG,iElem)/m_lMicro(iElem);
                   c_Tracc = c_Tracc+ ( m_NormalesMicro(:,:,iElem)'* ...
                       m_Sigma(:,iPG,iElem)/m_lMicro(iElem))*m_pesoPG(iPG);
                   
                   
                   pp = m_NormalesMicro(:,:,iElem)'*m_CT(:,:,iPG,iElem)*...
                    ( m_IndTens+m_BT(:,:,iPG,iElem)*m_VarFlucElem )/ m_lMicro(iElem);
                    
                   m_CTHomog_bb  = m_CTHomog_bb +...
                    (pp *...
                     m_NormalesMicro(:,:,iElem) / facNormMicro/m_lMicro(iElem) )* m_pesoPG(iPG);

                   m_CTHomog_ub  = m_CTHomog_ub + pp * m_pesoPG(iPG);
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
                end
             end
         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               c_Tracc{iSet} = m_Tracc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
      
%%%%%%%%%%%%%%%%%%%%%%%%%      
      m_SigmaHomog   = c_Tracc  /(facNormMicro*longFis);
      m_CTHomog_bb = m_CTHomog_bb /(facNormMicro*longFis);
      m_CTHomog_ub = m_CTHomog_ub /(facNormMicro*longFis);
      
      m_CTHomog(1:nDime,1:nDime)   = m_CTHomog_bb;
      m_CTHomog(nDime+1:2*nDime,:) = m_CTHomog_ub;
%%%%%%%%%%%%%%%%%%%%%%%%%      
%      m_SigmaHomog = f_HomogArea(c_Tracc,nDime,facNormMicro*longFis,c_DetJTLoc,e_DatSet,e_VG);
%      m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicroL,...
%          m_ElemLoc,m_ElemLoc,e_VG);
%%%%%%%%%%%%%%%%%%%%%%%%%      


      %Cuando el modelo constitutivo multiescala bifurca ya no se calcula m�s tensor tangente impl�cito, pero
      %igual se debe devolver una estructura, con cualquier valor para que el ensamblaje sea correcto en las
      %matrices m_TensorTang (ver como plantear esto de otra forma). No se puede descartar la parte que guarda
      %los tensores impl�citos ya que puede haber PG que ya bifurcaron y otros que no.
      if esImplexMacro
         m_CTHomog = struct('Implex',m_CTHomog,'Impli',zeros(ntens,ntens));
      end
   else
      %Se asume que previo a la bifurcaci�n y al imponer la discontinuidad fuerte, las tensiones old
      %homogeneizadas son iguales a las que se obtiene con la homogeneizaci�n de las sigma_old micro
      %de la celda unitaria (esto se realiza para ahorrarse la resta de las todas las tensiones
      %micro y homogeneizar el delta sigma).
%       m_DsigmaHomog = m_sigmaHomog-m_sigmaHomog_old;
  %    c_Dsigma = cellfun(@(x,y)minus(x,y),c_Tens_new,c_Tens_old,'UniformOutput',false);
      m_SigmaHomog = f_HomogArea(c_Tens_new,ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
      % Tensor Tangente Homogeneizado
%       m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,e_VG);  
      m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
         true(nElem,1),true(nElem,1),e_VG);
      %IMPLEX: Solo interesa que homogenice los tensores tangentes constitutivos impl�citos hasta 
      %la bifurcaci�n. Luego no se utiliza m�s, ya que la mismas se utilizan en el an�lisis de bifurcaci�n.
      if esImplexMacro
         %Se convierte c_CT con las matrices impl�citas micro �nicamente. Ver si no poner condicional implex
         %dentro de la funci�n f_ModTangHomog en su lugar (tambi�n hab�a que considerarla dentro de la funci�n
         %f_MatGlobal.
         for iSet = 1:nSet
            e_DatMatSet = e_DatSet(iSet).e_DatMat;
            nPG = e_DatSet(iSet).e_DatElem.npg;
            esImplex = e_DatMatSet.esImplex;
            if esImplex
               c_CT{iSet} = c_CT{iSet}(:,:,nPG+1:2*nPG,:);
            end            
         end
         %C�lculo de la matriz impl�cita
         KT = f_MatGlobal(c_CT,e_VarAux,e_DatSet,e_VG);
         %Homogenizaci�n del tensor impl�cito
         m_CTHomogImpli = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
            true(nElem,1),true(nElem,1),e_VG);
         m_CTHomog = struct('Implex',m_CTHomog,'Impli',m_CTHomogImpli);
      end     
   end
   
   %Se modifica tensor constitutivo tangente homogeneizado con el factor de proyecci�n que se utiliza para
   %inyectar las deformaciones macro.
   
   hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
      'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,'m_ElemLoc',m_ElemLoc,...
      'c_DefMacro',{c_DefMacro},'omegaMicroL',omegaMicroL,'lMacro',lMacro,'lMicro',lMicro,...
      'c_NormalesMicro',{c_NormalesMicro},'longFis',longFis,'facNormMicro',facNormMicro);     %,'m_TensProy',m_TensProy
   
end

function KT = f_MatGlobal(c_CT,e_VarAux,e_DatSet,e_VG)

      %Esta funci�n devuelve la matriz de rigidez global si se conoce cu�les son los tensores constitutivos
      %tangentes de todos los puntos de gauss de la estructura.
      nSet = e_VG.nSet;
      %ndn = e_VG.ndn;
      %
      c_Ke = cell(nSet,1);
      c_Fil = cell(nSet,1);
      c_Col = cell(nSet,1);
      %
      for iSet = 1:nSet
         
         %Ver si no es peor usar esto, por ejemplo que haga alguna copia innecesaria, que usar siempre
         %e_DatSet(iSet).
         e_DatiSet = e_DatSet(iSet);
         e_DatElemSet = e_DatiSet.e_DatElem;
         nPG = e_DatElemSet.npg;
         wg = e_DatElemSet.wg;
         dofpe = e_DatElemSet.dofpe;
         eltype = e_DatElemSet.eltype;
         nElem = e_DatiSet.nElem;
         %conec = e_DatiSet.conec;
         m_DofElem = e_DatiSet.m_DofElem;
         m_BT = e_DatiSet.m_BT;
         m_DetJT = e_DatiSet.m_DetJT;
         m_CT = c_CT{iSet};
         %
         m_Ke = zeros(dofpe*dofpe,nElem);     
         %
         %dofElemSet = f_DofElem(reshape(conec',1,[]),ndn);
         dofElemSet = m_DofElem(:);
         m_Fil = reshape(repmat(reshape(dofElemSet,dofpe,[]),dofpe,1),1,[]);
         m_Col = reshape(repmat(dofElemSet',dofpe,1),1,[]); 
         %
         switch eltype
            case {2,4,8,20}
               %Se est� considerando que en los elementos cuadrangulares 20 y 23, en los PGs arriba del 4 
               %vienen con pesos impuestos nulos. Por lo tanto �nicamente se hace la integraci�n de los 
               %valores que est�n sobre esos 4PGs.
               m_pesoPG = bsxfun(@times,m_DetJT,wg);
               %Con las matrices de deformaci�n precalculadas la determinaci�n de las matrices elementales
               %son iguales para los elementos est�ndars.
               for iElem = 1:nElem
               %parfor iElem = 1:nElem
                  m_kt = zeros(dofpe,dofpe);
                  for iPG = 1:nPG
                     B = m_BT(:,:,iPG,iElem);
                     ct = m_CT(:,:,iPG,iElem);
                     m_kt = m_kt+B'*ct*B*m_pesoPG(iPG,iElem);
                  end
                  %Para acelerar un poco el c�lculo se guarda en forma de columna las matrices de rigidez
                  %elementales.
                  m_Ke(:,iElem) = m_kt(:);
               end   
%             case 20
%                %El t�rmino constante de este elemento, proveniente del problema mixto, se integra con 1 solo
%                %punto de gauss, en el lugar de 4. Por ello al PG 5 se incorpora el peso 4, que corresponde al
%                %caso de un solo punto de gauss (wg viene con cero para realizar correctamente otras
%                %integraciones dentro del c�digo).
%                %wg(5) = 4;
%                m_pesoPG = bsxfun(@times,m_DetJT,wg);
%                % Variable auxiliar del elemento
%                %m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
%                for iElem = 1:nElem                  
%                %parfor iElem = 1:nElem
%                   %% VERSI�N 1                  
%                   %Condici�n de bifurcaci�n
%                   %Se determina mediante el an�lisis de bifurcaci�n del PG central (5), una vez que este
%                   %bifurc� todo el elemento se considera bifurcado (y en consecuencia los restantes PGs).
%                   condBif = m_VarAuxElem(1,iElem);                  
%                   % Par�metros de estabilizaci�n por PG
%                   %Se multiplica los par�metros de estabilizac�n por los pesos de los puntos de gauss directamente, as� tener
%                   %en cuenta este par�metro en la integraci�n. 
%                   if condBif>0
%                      %Si condBif es igual a 1, significa que se detect� bifurcaci�n.
%                      %Ver si no vale la pena activar la subintegraci�n unos pasos despu�s que bifurc� (condBif=2), o eso se
%                      %deja para la SD.
%                      m_pesoPGElem(1:4) = e_DatElemSet.estabBif*m_pesoPG(1:4,iElem);
%                      %Se asume que el PG central, donde se subintegra, est� en la posici�n 5 de la lista de puntos de gauss.
%                      m_pesoPGElem(5) = (1-e_DatElemSet.estabBif)*m_pesoPG(5,iElem); 
%                   elseif condBif==0
%                      %Si condBif es cero, significa que todav�a no se detect� la bifurcaci�n.
%                      m_pesoPGElem(1:4) = e_DatElemSet.estabNoBif*m_pesoPG(1:4,iElem);
%                      %Se asume que el PG central, donde se subintegra, est� en la posici�n 5 de la lista de puntos de gauss.
%                      m_pesoPGElem(5) = (1-e_DatElemSet.estabNoBif)*m_pesoPG(5,iElem); 
%                   end  
%                   m_kt = zeros(dofpe,dofpe);
%                   for iPG = 1:nPG
%                      B = m_BT(:,:,iPG,iElem);
%                      ct = m_CT(:,:,iPG,iElem);
%                      m_kt = m_kt+B'*ct*B*m_pesoPGElem(iPG);
%                   end
%                   %Para acelerar un poco el c�lculo se guarda en forma de columna las matrices de rigidez
%                   %elementales.
%                   m_Ke(:,iElem) = m_kt(:);
%                   %% VERSI�N 2
%                   %Condici�n de bifurcaci�n
%                   %Se determina mediante el an�lisis de bifurcaci�n del PG central (5), una vez que este
%                   %bifurc� todo el elemento se considera bifurcado (y en consecuencia los restantes PGs).
%                   condBif = m_VarAuxElem(1,iElem);    
%                   % Par�metros de estabilizaci�n por PG
%                   %Se multiplica los par�metros de estabilizac�n por los pesos de los puntos de gauss
%                   %directamente, as� tener en cuenta este par�metro en la integraci�n.
%                   if condBif>0
%                      %Si condBif es igual a 1, significa que se detect� bifurcaci�n.
%                      %Ver si no vale la pena activar la subintegraci�n unos pasos despu�s que bifurc� (condBif=2), o eso se
%                      %deja para la SD.
%                      paramEstab = e_DatElemSet.estabBif;
%                   elseif condBif==0
%                      %Si condBif es cero, significa que todav�a no se detect� la bifurcaci�n, o que descarg�.
%                      paramEstab = e_DatElemSet.estabNoBif;
%                   end
%                   paramEstabRest = 1-paramEstab;
%                   B5 = m_BT(:,:,5,iElem); 
%                   ct5Estab = paramEstabRest*m_CT(:,:,5,iElem);
%                   m_kt = zeros(dofpe,dofpe);
%                   for iPG = 1:nPG-1
%                      B = m_BT(:,:,iPG,iElem);
%                      ctEstab = paramEstab*m_CT(:,:,iPG,iElem);
%                      m_kt = m_kt+(B'*ctEstab*B+B5'*ct5Estab*B5)*m_pesoPG(iPG,iElem);
%                   end
%                   %Para acelerar un poco el c�lculo se guarda en forma de columna las matrices de rigidez
%                   %elementales.
%                   m_Ke(:,iElem) = m_kt(:);
%                   %% VERSI�N 3
%                   %Se asume que las matrices m_CT son las impl�citas estabilizadas, por lo que se integra en
%                   %el elemento con cuatro punto de gauss (se podr�a integrar en los 5 puntos de gauss ya que
%                   %el punto de gauss 5 tiene peso nulo).
%                   m_kt = zeros(dofpe,dofpe);
%                   for iPG = 1:nPG-1
%                      B = m_BT(:,:,iPG,iElem);
%                      ct = m_CT(:,:,iPG,iElem);
%                      m_kt = m_kt+B'*ct*B*m_pesoPG(iPG,iElem);
%                   end
%                   %Para acelerar un poco el c�lculo se guarda en forma de columna las matrices de rigidez
%                   %elementales.
%                   m_Ke(:,iElem) = m_kt(:);
%                end   
            otherwise
               error(['Modelo Multiescala Cohesivo: Matriz de rigidez global impl�cita para el implex: ',...
                  'Tipo de elemento no definido'])
         end
         c_CT{iSet} = m_CT;
         c_Ke{iSet} = m_Ke(:);
         c_Fil{iSet} = m_Fil;
         c_Col{iSet} = m_Col;
         
      end
      
      % Ensamble de matriz de rigidez global
      KT = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:})); 
      
end