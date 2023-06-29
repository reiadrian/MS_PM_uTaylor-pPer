function [m_CTHomog,sigmaHomog,hvar_newMacro] = f_RMap_ME(...
   eps_new,hvar_oldMacro,e_DatMatSetMacro,e_VGMacro)

   %Se recupera variables micro
   xx = e_DatMatSetMacro.xx;
   omegaMicro = e_DatMatSetMacro.omegaMicro;
   e_DatSet = e_DatMatSetMacro.e_DatSet;
   e_VG = e_DatMatSetMacro.e_VG;
   esImplexMacro = e_DatMatSetMacro.esImplex;
   e_VarEst_old = hvar_oldMacro.e_VarEst;
   u = hvar_oldMacro.u;
   c_GdlCond = hvar_oldMacro.c_GdlCond;
   Fint         = hvar_oldMacro.Fint;
   m_LinCond    = hvar_oldMacro.m_LinCond;
   doff         = hvar_oldMacro.doff;
   dofl         = hvar_oldMacro.dofl;
   %c_DefMacro   = hvar_oldMacro.c_DefMacro;

   %vfix = hvar_oldMacro.vfix;
   e_VarAux = hvar_oldMacro.e_VarAux;

   % VARIABLES GLOBALES
   nElemTot = e_VG.nElem;
   ntens = e_VG.ntens;
   ndoft = e_VG.ndoft;

   % INICIALIZACION DE VARIABLES
   Du_step_old = zeros(ndoft,1);  % Vector de incrementos de desplazamientos en el paso de tiempo previo
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
   
    e_VG.SharedFS = e_VGMacro.SharedFS; %AA
   
   %Nombre interno de la celda unitaria 
   %(hacer esto en cada iteraci�n puede ser medio lento, ver que hacer sino).
   %e_VG.fileCompleto = [e_VG.fileCompleto,'_EM',int2str(e_VGMacro.iElemNum),'PGM',int2str(e_VGMacro.iPG)];
   %Se guarda los datos del matlabpool (se utiliza para imprimir de pasos y iteraciones a nivel micro que no
   %convergieron).
   e_VG.nLab = e_VGMacro.nLab; %AA
   e_VG.tipoMPool = e_VGMacro.tipoMPool; %AA

   % DESPLAZAMIENTO IMPUESTO
   %No ser�a necesario estas operaciones porque vfix en todos los modelos cl�sicos de las
   %formulaciones multiescala son nulos. Adem�s que habr�a que interpretar como realizar el delta
   %psi_value a nivel micro, ya deber�a se corresponder con el incremento de tiempo a nivel macro,
   %pero lo que implicar�a que en cada paso de tiempo se est� resolviendo un RVE distinto.
   %vDeltaFixTemp = vfix;
   %vDeltaFixTemp(doffCondCte) = vfix(doffCondCte)*(psi_value - psi_value_old);
   %u(doff) = u(doff) + m_InvCRR*vDeltaFixTemp(doff);
  
   % Deformaci�n macro aplicada en la Celda unitaria.
   %Se aplica en forma uniforme en todo dominio.
   %Por simplificidad se considera una deformaci�n macro aplicada distinta por elemento, y no por PG.
   %(no es necesario considerar una estructura para esta variable).
   %m_DefMacro = repmat(eps_new,1,nElem);
   %Deformaci�n macro por elemento y por punto de gauss, dividida en sets.
   c_DefMacro = arrayfun(@(x)repmat(eps_new,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false); % AA:Desnloquee
   %Deformaci�n macro por elemento, dividad en sets.
%    c_DefMacro = arrayfun(@(x)repmat(eps_new,[1,x.nElem]),e_DatSet,'UniformOutput',false);
   % AA: Bloeque� linea 68 c_DefMacro = arrayfun(@(x)repmat(eps_new,[1,x.nElem]) 
   %ticIDNewt = tic;
   % ESQUEMA DE NEWTON-RAPHSON
   Du_step_new = zeros(ndoft,1);
   [u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] = newton_raphson(...
      xx,m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
      c_DefMacro,e_VG);
  
   %Se guarda las fuerzas internas s�lo para impresi�n de las gr�ficas X-Y
   %Se est� usando equilibrio en totales, por lo que no es necesario la siguiente l�nea.
   %Fint = Fint+DFint;
   %
   %fprintf('Tiempo del Newton micro: %f\n',toc(ticIDNewt));   
   
   % TENSOR TANGENTE HOMOGENEIZADO
   m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
        true(nElemTot,1),true(nElemTot,1),e_VG);
   
   %Se asume que no se realiza an�lisis de bifurcaci�n con el tensor tangente constitutivo homogeneizado, por
   %lo que en el caso ser implex, se devuelve nulo el tensor impl�cito homogeneizado.
   if esImplexMacro
      m_CTHomog = struct('Implex',m_CTHomog,'Impli',zeros(ntens,ntens));
   end
   
   % C�LCULO DE VARIABLES HOMOGENEIZADAS
   sigmaHomog = f_HomogArea({e_VarEst_new.sigma},ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
   %defHomog = f_HomogArea({e_VarEst_new.eps},ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
   defHomogFl = f_HomogArea({e_VarEst_new.eps_fluct},ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
   %Verificaci�n de la media del desplazamiento
   m_uMedioFl = f_MediaDespCU(u,omegaMicro,e_DatSet,e_VG);
   %
   fprintf('Elemento %d: PG %d: Norma de la deformaci�n fluctuante media: %g\n',e_VGMacro.iElemNum,...
      e_VGMacro.iPG,norm(defHomogFl))
   fprintf('Elemento %d: PG %d: Norma del desplazamiento fluctuante medio: %g\n',e_VGMacro.iElemNum,...
      e_VGMacro.iPG,norm(m_uMedioFl))
   
   hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
      'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,'c_DefMacro',{c_DefMacro});
   
end
