function hVarNew = f_DefinicionhVar(conshyp,sihvarpg,nPG,protype)

   %Esta funcion inicializa las variables historicas segun el modelo constitutivo, considerando que el caso
   %multiescala puede ser una array de estructura y para el estandar una matriz (es para llamarse dentro de la
   %funcion del elemento). Esta generacion evita transferencia de datos a los nodos (habria que ver el tiempo
   %adicional agregado por el llamado de esta funcion).
   %sihvarpg = e_DatMatSet.sihvarpg;
   %nPG = e_DatElemSet.npg;
   %conshyp = e_DatMatSet.conshyp;   
   switch conshyp
       %#######################################################################
       case {1,2,4,5,8,10,11,12,13,14,15,16,17,52,100,110} % AA: add case 17
         hVarNew = zeros(sihvarpg,nPG);
         %#######################################################################
       case {50,55,60,61,62}
           if protype==0
               hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
                   'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);
           elseif protype==1
                %###################################################################################
               %Agregue Fext
               hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
                   'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[],'c_GradPorMacro_dup',[],'c_GradPorMacro_up_n',[],...
                   'c_PorMacro',[],'c_DefMacro_new',[],'c_GradPorMacro_up_new',[],...
                   'c_PorMacro_new',[],'Fext',[]);              
               %###################################################################################
           end
       case 51
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
            'm_LinCond',[],'doff',[],'dofl',[],'m_ElemLoc',[],'c_DefMacro',[],'omegaMicroL',[],...
            'lMacro',[],'lMicro',[],'c_NormalesMicro',[],'longFis',[],'facNormMicro',[]);        %,'m_TensProy',[]
       case 53
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
            'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]
       case 54
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],...
             'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],...
             'doff',[],'dofl',[],'m_ElemLoc',[],'c_DefMacro',[],...
             'omegaMicroL',[],'lMacro',[],'lMicro',[],'c_NormalesMicro',[],...
             'longFis',[],'facNormMicro',[] );        %,'m_TensProy',[]        
        otherwise
         error('Matrices Elementales: Variables Historicas: Inicializacion: Modelo constitutivo no definido.')         
   end
   
end