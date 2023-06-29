function f_PostProcME_Bif(iStep,m_NumElemMacro,m_DefMacro,m_phiMacro,m_porpMacro,...
    e_DatMatSetMacro,hvar_newMacro,e_VGMacro)


   %Se recupera variables micro
   xx = e_DatMatSetMacro.xx;
   in = e_DatMatSetMacro.in;
   e_VG = e_DatMatSetMacro.e_VG;
   m_ElemPGImpr = e_DatMatSetMacro.m_ElemPGImpr;
   e_DatSet = e_DatMatSetMacro.e_DatSet;
   m_SetElem = e_DatMatSetMacro.m_SetElem;
   %Se guarda el paso tiempo a nivel micro.
   e_VG.istep = e_VGMacro.istep;
   %
   postpro_impre_step = e_VG.postpro_impre_step;
   IRES = e_VG.IRES;
   ndime = e_VG.ndime;
   ndn_d = e_VG.ndn_d;
   ndn_p = e_VG.ndn_p;
   nElem = e_VG.nElem;
   %Nombre completo del archivo micro
   nomArchCompl = e_VG.fileCompleto;
   %   
   for iPGGraf = 1:size(m_ElemPGImpr,2)
      
      iElemSetMacro = m_ElemPGImpr(1,iPGGraf);
      iPGMacro = m_ElemPGImpr(2,iPGGraf);
      %
      %m_iDefMacro = m_DefMacro(:,iPGMacro,iElemSetMacro);
      m_iDefMacro = m_DefMacro(:,iElemSetMacro);
      m_iphiMacro = m_phiMacro(:,iElemSetMacro);
      m_iporpMacro = m_porpMacro(:,iElemSetMacro);
      e_VarEst_new = hvar_newMacro(iPGMacro,iElemSetMacro).e_VarEst;
      e_VarAux = hvar_newMacro(iPGMacro,iElemSetMacro).e_VarAux;
      u = hvar_newMacro(iPGMacro,iElemSetMacro).u;
      c_GdlCond = hvar_newMacro(iPGMacro,iElemSetMacro).c_GdlCond;
      Fext = hvar_newMacro(iPGMacro,iElemSetMacro).Fext;
      Fint = hvar_newMacro(iPGMacro,iElemSetMacro).Fint;
      %Deformaciones MACRO-APPLIED al paso de tiempo "n+1"
      m_DefMacroCU = hvar_newMacro(iPGMacro,iElemSetMacro).c_DefMacro;
      %Gradiente de poro presiones MACRO-APPLIED al paso de tiempo "n+1"
      m_GradPorMacroCU = hvar_newMacro(iPGMacro,iElemSetMacro).c_GradPorMacro_dup;
      %
      m_ElemLoc = false(1,nElem);
      %Se cambia el nombre del archivo de postproceso de GiD de la celda unitaria, para
      %individualizar el elemento y el punto de Gauss que corresponde.
      e_VG.fileCompleto = [nomArchCompl,'_EM',num2str(m_NumElemMacro(iElemSetMacro)),'PGM',...
         num2str(iPGMacro)];

      % IMPRESION DE RESULTADOS
      %Notar que se imprime este paso si a nivel macro este paso se imprime, y a la vez, que si a
      %nivel micro debe imprimirse.
      if IRES   
         if ~mod(iStep,postpro_impre_step)
            % DESPLAZAMIENTO TOTAL DE LA MICRO-CELDA
            ud=u(e_VG.pos_dG);
            udTotal = [m_iDefMacro(1),m_iDefMacro(4)/2;m_iDefMacro(4)/2,m_iDefMacro(2)]*xx(:,1:2)'...
               +reshape(ud,ndn_d,[]);
           %############################################################################################
            u = f_porp_nint(u,e_DatSet,e_VG); %AA: obtiene los valores de poropresiones en los nodos internos
            %############################################################################################
            up=u(e_VG.pos_pG);
            ndn_p=e_VG.ndn_p;
            upTotal = m_iporpMacro+m_iphiMacro'*xx(:,1:2)'+reshape(up,ndn_p,[]);
            % POSTPROCESO GiD
            matlab2gid_res_Bif(iStep,in,u,c_GdlCond,e_DatSet,e_VarEst_new,e_VarAux,...
                 m_DefMacroCU,m_GradPorMacroCU,udTotal,upTotal,Fext,Fint,m_ElemLoc,e_VG)
         end
      end
      %############################################################################################
      if e_VG.protype==1 %AA
         u(3*e_VG.in_int) = 0.0; %AA: anulo las poropresiones en los lados (nodos internos)
      end %protype  
      %############################################################################################
      %Impresion graficas
      sigmaHomog = e_VGMacro.sigmaHomog; %AA
      epsilon_Macro = e_VGMacro.epsilon_Macro; %AA: agregue ambas variables para poder ingresarlas con valores
      f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,e_DatSet,m_SetElem,sigmaHomog,epsilon_Macro,e_VG) %AA
      
   end
         
end