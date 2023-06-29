function e_VarEst = f_eVarEstInic(c_Field,e_DatSet,e_VG,varargin)
   
%Inicializa la estructuras de variable de estado
   nSet = e_VG.nSet;
   ntens = e_VG.ntens;
%    ndoft = e_DatSet.e_DatMat.e_VG.ndoft;
   ndime = e_VG.ndime;
   nField = length(c_Field);
   c_StructDef = cat(1,c_Field,repmat({[]},1,nField));
   %e_VarEst(1:nSet,1) = struct(c_StructDef{:});
   %Si son estructuras vacias tambien se puede usar
   e_VarEst(nSet,1) = struct(c_StructDef{:});
   for iSet = 1:nSet
      sitvare = e_DatSet(iSet).sitvare;
      % AA22: Recupero informacion de Nro de PG para usar velflu y mflu 
      % linea 190 en adelante
      npg = sitvare/ntens; 
      sihvare = e_DatSet(iSet).sihvare;
      nElem = e_DatSet(iSet).nElem;
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      conshyp = e_DatMatSet.conshyp;
      nVarHistElem = e_DatElemSet.nVarHistElem;

      m_InicZerosSet = zeros(sitvare,nElem);
      for iField = 1:nField
         if strncmp(c_Field{iField},'hvar',4)
            switch conshyp
                %####################################################################################
               case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17.18,52,100,110} %AA: case 17
                   %####################################################################################
                  v_InicHvarSet = zeros(sihvare,nElem);
                  %Se inicializa la variable historica para el modelo constitutivo 11
                  if conshyp==11
                     %Se inicializa la variables r_old(5), q_old(6), rInic_old(8) y qInic_old(9). Todas tienen
                     %que ser igual a r0.
                     sihvarpg = e_DatMatSet.sihvarpg;
                     if e_DatMatSet.tit==1
                        %En el caso de ser exponencial tambien se inicializa las rInic_old(end-1) y
                        %qInic_old(end).                        
                        m_Ind = [5;6;sihvarpg-1;sihvarpg];
                     else
                        m_Ind = [5;6];
                     end
                     v_InicHvarSet(bsxfun(@plus,m_Ind,0:sihvarpg:sihvare-1),:) = e_DatMatSet.r_0;                       
                  elseif conshyp==110
                     sihvarpg = e_DatMatSet.sihvarpg;
                     m_Ind = [3;4;5];
                     v_InicHvarSet(bsxfun(@plus,m_Ind,0:sihvarpg:sihvare-1),:) = 1;                       
                  end                      
               case {50,55,60,61,62}   %Modelo MultiEscala continuo
                  clear v_InicHvarSet
                  if nargin>3&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el calculo de las matrices de
                     %condiciones de borde, que principalmente las periodicas son muy lentas, y no se
                     %utilizan.
                     if e_VG.protype == 0
                         v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);
                     elseif e_VG.protype == 1
                         %###################################################################################
                         %Agregue Fext
                         v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[],...
                        'c_GradPorMacro_dup',[],'c_GradPorMacro_up_n',[],'c_PorMacro',[],'Fext',[],...
                        'c_DefMacro_new',[],'c_GradPorMacro_up_new',[],'c_PorMacro_new',[]);
                        %###################################################################################
                     end
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     %Como es el mismo set, se asume que las celdas unitarias de todos los puntos
                     %de gauss de todos los elementos tiene las mismas condiciones de borde iniciales.
                     %Como las condiciones de borde luego puede variar segun la celda (si bifurca o
                     %no), se las almacena como variables historicas (en lugar en las e_DatSet, que se
                     %asume que son datos fijos).
                     %Se descarta todo las matrices vfix, m_InvCRR y doffCondCte que se refieren a los
                     %desplazamientos impuestos de grados de libertad fijos, porque se asumen que son
                     %nulos. Ver si hay alguna tipo de celda unitaria donde se necesario considerar
                     %algo distinto a esto.
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                        e_DatMatSet.xx,e_DatSetMicro,e_VGMicro.m_ConecFront);
                     %Las variables u (desplazamiento fluctuante micro), c_GdlCond (variable condensada
                     %micro) y e_VarAux (variables auxiliares micro) se hacen variables tipo historicas
                     %ya que en la celda unitaria se aplica una deformacion macro variable en cada
                     %iteracion macro mientras que en la iteracion micro se consigue que siempre se
                     %parta de la misma condicion inicial u (fluctuante micro).
                     %Otra posibilidad es que esta variables se vayan actualizando en cada iteracion
                     %macro, es decir pasandola como variable auxiliar macro (tambien se podria
                     %conseguir el mismo efecto si se pasa como argumento las variables new macro, y se
                     %utiliza los valores obtenidos de la misma). Esto haria que una mala iteracion
                     %macro, que se aleje de la solucion, haga que en la proxima iteracion macro se
                     %parta en la iteracion micro de una condicion inicial capaz muy alejada de la
                     %solucion (aunque en forma opuesta, si las iteraciones macro son correctas, si se
                     %utiliza condiciones iniciales actualizadas, se necesitaria capaz menos
                     %iteraciones micro). Ahora almacenandola como variable historica y usando los
                     %valores obtenidos del old macro, siempre se parte de una condicion inicial micro
                     %que se sabe convergido en el paso previo, pero que puede exigir mas iteraciones
                     %micro ya que puede estar mas alejado de la solucion.
                     if e_VG.protype == 0
                         v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatMatSet.xx,e_DatSetMicro,e_VGMicro),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.nElem),e_DatSetMicro,...
                           'UniformOutput',false)});
                     elseif e_VG.protype == 1
                         v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatMatSet.xx,e_DatSetMicro,e_VGMicro),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.nElem),e_DatSetMicro,...
                           'UniformOutput',false)},...
                           'c_GradPorMacro_dup',{arrayfun(@(x)zeros(2,x.nElem),e_DatSetMicro,'UniformOutput',false)},...
                           'c_GradPorMacro_up_n',{arrayfun(@(x)zeros(2,x.nElem),e_DatSetMicro,'UniformOutput',false)},...                           
                           'c_PorMacro',{arrayfun(@(x)zeros(1,x.nElem),e_DatSetMicro,'UniformOutput',false)},...
                           'c_DefMacro_new',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.nElem),e_DatSetMicro,...
                           'UniformOutput',false)},...
                           'c_GradPorMacro_up_new',{arrayfun(@(x)zeros(2,x.nElem),e_DatSetMicro,'UniformOutput',false)},...                           
                           'c_PorMacro_new',{arrayfun(@(x)zeros(1,x.nElem),e_DatSetMicro,'UniformOutput',false)},...
                           'm_VarFluc_eps0',zeros(e_DatSet.e_DatMat.e_VG.ndoft,ntens),...
                           'm_VarFluc_p0',zeros(e_DatSet.e_DatMat.e_VG.ndoft,1),...
                           'm_VarFluc_phi0',zeros(e_DatSet.e_DatMat.e_VG.ndoft,2));
                      end
                        %Notar que la deformacion para grandes deformaciones (Gradiente de deformacion F) no
                        %se inicializa en la identidad ya que este es pisado siempre en la funcion f_RMap_ME.
                        %{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),e_DatSetMicro,'UniformOutput',false)});
                  end
               case 51   %Modelo MultiEscala con discontinuidad fuerte (fisura cohesiva)
                  clear v_InicHvarSet
                  if nargin>3&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el calculo de las matrices de
                     %condiciones de borde, que principalmente las periodicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'m_ElemLoc',[],...
                        'c_DefMacro',[],'omegaMicroL',[],'lMacro',[],'lMicro',[],...
                        'c_NormalesMicro',[],'longFis',[],'facNormMicro',[]);    %,'m_TensProy',[]
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                        e_DatMatSet.xx,...
                        e_DatSetMicro,e_VGMicro.m_ConecFront);
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatMatSet.xx,e_DatSetMicro,e_VGMicro),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'm_ElemLoc',[],...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                           e_DatSetMicro,'UniformOutput',false)},...
                        'omegaMicroL',e_DatMatSet.omegaMicro,'lMacro',0,'lMicro',0,...
                        'c_NormalesMicro',{cell(e_VGMicro.nSet,2)},'longFis',0,'facNormMicro',1);   %'m_TensProy',[]
                  end
               case 53   %Modelo MultiEscala BCNA, SANTA FE
                  clear v_InicHvarSet
                  if nargin>3&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el calculo de las matrices de
                     %condiciones de borde, que principalmente las periodicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                        'c_DefMacro',[]);    %,'m_TensProy',[]
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                        e_DatMatSet.xx,...
                        e_DatSetMicro,e_VGMicro.m_ConecFront);
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                           e_DatSetMicro,'UniformOutput',false)});   %'m_TensProy',[]
                  end
               case 54   %Modelo MultiEscala BCNA, SANTA FE
                  clear v_InicHvarSet
                  if nargin>3&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el calculo de las matrices de
                     %condiciones de borde, que principalmente las periodicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                        'm_ElemLoc',[],...
                        'c_DefMacro',[],'omegaMicroL',[],'lMacro',[],'lMicro',[],...
                        'c_NormalesMicro',[] ,...
                        'longFis',[],'facNormMicro',[] );   %'m_TensProy',[]
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                        e_DatMatSet.xx,...
                        e_DatSetMicro,e_VGMicro.m_ConecFront);
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'm_ElemLoc',[],...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                           e_DatSetMicro,'UniformOutput',false)}, ...   %'m_TensProy',[]
                        'omegaMicroL',e_DatMatSet.omegaMicro,'lMacro',0,'lMicro',0,...
                        'c_NormalesMicro',{cell(e_VGMicro.nSet,1)},...
                        'longFis',0,'facNormMicro',1 );
                  end
               otherwise
                  error('Inicializacion variables historicas: Modelo constitutivo no definido.')
            end
            e_VarEst(iSet).(c_Field{iField}) = v_InicHvarSet;
         else
            if strncmp(c_Field{iField},'phi',3)
                e_VarEst(iSet).(c_Field{iField}) = zeros(ndime*npg,nElem);
%             elseif strncmp(c_Field{iField},'porpr',5)
%                 e_VarEst(iSet).(c_Field{iField}) = zeros(4,nElem);
                %Representa las poropresiones nodales
            elseif strncmp(c_Field{iField},'velflu',6)
                e_VarEst(iSet).(c_Field{iField}) = zeros(ndime*npg,nElem);
            elseif strncmp(c_Field{iField},'mflu',4)
                e_VarEst(iSet).(c_Field{iField}) = zeros(npg,nElem);
            elseif strncmp(c_Field{iField},'velflu_sta',10)
                e_VarEst(iSet).(c_Field{iField}) = zeros(ndime*npg,nElem);
            elseif strncmp(c_Field{iField},'velflu_total',12)
                e_VarEst(iSet).(c_Field{iField}) = zeros(ndime*npg,nElem);
            elseif strncmp(c_Field{iField},'phi_fluct',9)
                e_VarEst(iSet).(c_Field{iField}) = zeros(ndime*npg,nElem);  
            elseif strncmp(c_Field{iField},'p_fluct',7)
                e_VarEst(iSet).(c_Field{iField}) = zeros(4,nElem);
            elseif strncmp(c_Field{iField},'porpr',5)
                e_VarEst(iSet).(c_Field{iField}) = zeros(npg,nElem);
                %Representa las poropresiones fluctuantes nodales
            else
                e_VarEst(iSet).(c_Field{iField}) = m_InicZerosSet;
            end
         end
      end
      %
      %Variables historica del elemento
      e_VarEst(iSet).VarHistElem = zeros(nVarHistElem,nElem);
   end
   
end