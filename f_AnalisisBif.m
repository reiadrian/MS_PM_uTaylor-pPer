function hvar_new = f_AnalisisBif(eps_new,hvar_new,e_DatSet,e_VG)
   
   ntens = e_VG.ntens;
   nElem = e_DatSet.nElem;
   e_DatElemSet = e_DatSet.e_DatElem;
   e_DatMatSet = e_DatSet.e_DatMat;
   nPG = e_DatElemSet.npg;
   m_Ce = e_DatMatSet.ce;
   poiss = e_DatMatSet.poiss;
   E = e_DatMatSet.young;
   gfv = e_DatMatSet.gfv;
   r_0 = e_DatMatSet.r_0;
   ftult = e_DatMatSet.ftult;
   q_0 = r_0;
   tita = e_DatMatSet.tit;
   esImplex = e_DatMatSet.esImplex;   
   %
   %Se regulariza el modelo constitutivo seg�n el tama�o de elemento.
   m_hReg = e_DatMatSet.m_hReg;
   switch tita
      case 0
         % Lineal            
         if e_DatMatSet.young~=0
            %H = -Sigmau^2/2/E/Gf*h;
            m_hba = -ftult^2/2/E/gfv*m_hReg;
         else
            m_hba = zeros(1,nElem);
         end
      case 1
         % Exponencial         
         %Para considerar que en el caso se utilice E=0, y el modelo constitutivo funcione
         %(igualmente hay que resolver que no queden nodos sin rigidez)
         if e_DatMatSet.young~=0
            %H = -Sigmau^2/E/Gf*h;
            m_hba = -ftult^2/E/gfv*m_hReg;
         else
            m_hba = zeros(1,nElem);
         end
      otherwise
         error('Modelo de da�o solo tracci�n regularizado: Modelo de ablandamiento no definido.')
   end
   
   eps_new = reshape(eps_new,ntens,nPG,nElem);
   m_SigmaEfect = zeros(ntens,nPG,nElem);
   for iPG = 1:nPG
      m_SigmaEfect(:,iPG,:) = m_Ce*reshape(eps_new(:,iPG,:),ntens,nElem);
   end
   % Determinaci�n de los autovalores positivos vectorizada para versi�n 2D
   %Para la determinaci�n de los valores principales se resuelve �nicamente en el plano xy (se ignora la
   %sigma_zz). Tambi�n que se tiene [sigma_xx,sigma_yy,sigma_zz,sigma_xy].
   if e_VG.ndime==2
      m_MediaTensNorm = sum(m_SigmaEfect(1:2,:,:),1)/2;
      m_DifTensNorm = m_SigmaEfect(1,:,:)-m_SigmaEfect(2,:,:);
      m_Raiz = hypot(m_DifTensNorm/2,m_SigmaEfect(4,:,:));
      m_TensPrinc1 = m_MediaTensNorm+m_Raiz;
      m_TensPrinc2 = m_MediaTensNorm-m_Raiz;
      m_AngDirP1 = atan2(2*m_SigmaEfect(4,:,:),m_DifTensNorm)/2;
      %Se toma la parte positiva del tensor de tensiones.
      %Ver si no hacer reshape para que sea una matriz 2D, ya que se est� trabajando con matrices
      %(1,nPG,nElem).
      m_TensPrinc1(m_TensPrinc1<0) = 0;
      m_TensPrinc2(m_TensPrinc2<0) = 0;
   else
      error('Da�o solo tracci�n: An�lisis de bifurcaci�n: No definido los autovalores positivo en el caso 3D.')
   end
   
   %Se asume un caso 2D.
   if e_VG.ndime==2
      hvar_new = reshape(hvar_new,[],nPG,nElem); 
      % Determinaci�n del criterio de bifurcaci�n      
      %Ver si hacer este an�lisis solo para los PG no bifurcados.
      m_r = hvar_new(1,:,:);
      m_q = hvar_new(2,:,:);
      m_DanioMas1 = m_q./m_r;
      m_ZCrit = (1+poiss)/E./m_DanioMas1.*((1-poiss)*(m_TensPrinc1+m_TensPrinc2).^2-2*m_TensPrinc1.*m_TensPrinc2);
      m_HCrit = m_DanioMas1-m_r.^2./m_ZCrit;
      if esImplex
         m_PGBifOld = hvar_new(6,:,:);
      else
         m_PGBifOld = hvar_new(5,:,:);
      end
      m_PGNoBifOld = ~m_PGBifOld;
      %
      if tita==0  %Lineal
         m_PGBif = bsxfun(@le,reshape(m_hba,1,1,nElem),m_HCrit);
      elseif tita==1  %Exponencial
         m_hba = repmat(reshape(m_hba,1,1,nElem),[1,nPG,1]);
         m_Hr = m_hba.*exp(m_hba.*(m_r-r_0)/q_0);
         m_PGBif = m_Hr<=m_HCrit;
      end
      %Por si acaso, si despu�s de la bifurcaci�n cambia la condici�n de bifurcaci�n, se mantiene en bifurcado
      %si ya se la detect� previamente.
      m_PGBif = m_PGBif|m_PGBifOld;
      %
      % Determinaci�n del �ngulo de bifurcaci�n
%       facPoiss = poiss/(1-poiss);
%       %Ver si directamente solo calcular las tensiones principales de los PG no bifurcados.
%       m_TangCuad = -(m_TensPrinc2(m_PGNoBifOld)-facPoiss*m_TensPrinc1(m_PGNoBifOld))./...
%          (m_TensPrinc1(m_PGNoBifOld)-facPoiss*m_TensPrinc2(m_PGNoBifOld));
%       %En el caso de ser negativo m_TangCuad los �ngulos cr�ticos son nulos (el �ngulo cr�tico es igual al
%       %�ngulo principal 1), y si el resultado es NaN (posiblemente por ser
%       %m_TensPrinc1-facPoiss*m_TensPrinc2=0 y m_TensPrinc2-facPoiss*m_TensPrinc1=0, 0/0 est� indefinido, que
%       %como m_TensPrinc1>m_TensPrinc2 y facPoiss>0, solo pasa cuando m_TensPrinc1=0 y m_TensPrinc2=0).
%       m_TangCuad(m_TangCuad<0|isnan(m_TangCuad)) = 0;
%       m_AngCrit = atan(sqrt(m_TangCuad));      
%       if esImplex
%          hvar_new(4,m_PGNoBifOld) = m_AngDirP1(m_PGNoBifOld)+m_AngCrit;
%          hvar_new(5,m_PGNoBifOld) = m_AngDirP1(m_PGNoBifOld)-m_AngCrit;
%          hvar_new(6,:,:) = m_PGBif;
%       else
%          hvar_new(3,m_PGNoBifOld) = m_AngDirP1(m_PGNoBifOld)+m_AngCrit;
%          hvar_new(4,m_PGNoBifOld) = m_AngDirP1(m_PGNoBifOld)-m_AngCrit;
%          hvar_new(5,:,:) = m_PGBif;
%       end
   %
      % Determinaci�n del �ngulo de bifurcaci�n
      %Ac� se est� recalculando los �ngulos de bifurcaci�n siempre aunque ya est� bifurcacado, esto puede que
      %este �ngulo cambie despu�s de la bifurcaci�n.
      facPoiss = poiss/(1-poiss);
      %Ver si directamente solo calcular las tensiones principales de los PG no bifurcados.
      m_TangCuad = -(m_TensPrinc2-facPoiss*m_TensPrinc1)./(m_TensPrinc1-facPoiss*m_TensPrinc2);
      %En el caso de ser negativo m_TangCuad los �ngulos cr�ticos son nulos (el �ngulo cr�tico es igual al
      %�ngulo principal 1), y si el resultado es NaN (posiblemente por ser
      %m_TensPrinc1-facPoiss*m_TensPrinc2=0 y m_TensPrinc2-facPoiss*m_TensPrinc1=0, 0/0 est� indefinido, que
      %como m_TensPrinc1>m_TensPrinc2 y facPoiss>0, solo pasa cuando m_TensPrinc1=0 y m_TensPrinc2=0).
      m_TangCuad(m_TangCuad<0|isnan(m_TangCuad)) = 0;
      m_AngCrit = atan(sqrt(m_TangCuad));      
      if esImplex
         hvar_new(4,:,:) = m_AngDirP1+m_AngCrit;
         hvar_new(5,:,:) = m_AngDirP1-m_AngCrit;
         hvar_new(6,:,:) = m_PGBif;
      else
         hvar_new(3,:,:) = m_AngDirP1+m_AngCrit;
         hvar_new(4,:,:) = m_AngDirP1-m_AngCrit;
         hvar_new(5,:,:) = m_PGBif;
      end 
   else
      error('Da�o solo tracci�n: An�lisis de bifurcaci�n: No definido el �ngulo de bifurcaci�n en el caso 3D.')
   end
   
   hvar_new = reshape(hvar_new,[],nElem);

end