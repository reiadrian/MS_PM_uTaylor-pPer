function [ct,sigma_new,eps_new,hvar_new,aux_var] = RMapDanoSTraccReg(...
   eps_new,hvar_old,aux_var,e_DatMatSet,e_VG)

   %Se regulariza el modelo constitutivo seg�n el tama�o de elemento.
   tita = e_DatMatSet.tit;
   %Ver si no pasar la matriz hReg de otra forma, ya que depende expl�citamente que se pase en el
   %e_VG el n�mero de elemento que le corresponde en el set, algo se hab�a tratado de evitar a este
   %nivel (solo se usaba para debug o imprensi�n).
   %Tambi�n se podr�a interpretar que es una propiedad que depende de la posici�n donde se eval�a, y
   %donde se asume que la posici�n se indica como el n�mero de elemento finito.
   hElemReg = e_DatMatSet.m_hReg(e_VG.iElemSet);
   switch tita
      case 0
         % Lineal            
         if e_DatMatSet.young~=0
            %H = -Sigmau^2/2/E/Gf*h;
            e_DatMatSet.hba = -e_DatMatSet.ftult^2/2/e_DatMatSet.young/e_DatMatSet.gfv*hElemReg;
         else
            e_DatMatSet.hba = 0;
         end
         [ct,sigma_new,eps_new,hvar_new,aux_var] = RMapDanoSoloTraccion(eps_new,hvar_old,aux_var,...
            e_DatMatSet,e_VG);
      case 1
         % Exponencial         
         %Para considerar que en el caso se utilice E=0, y el modelo constitutivo funcione
         %(igualmente hay que resolver que no queden nodos sin rigidez)
         if e_DatMatSet.young~=0
            %H = -Sigmau^2/E/Gf*h;
            e_DatMatSet.hba = -e_DatMatSet.ftult^2/e_DatMatSet.young/e_DatMatSet.gfv*hElemReg;
         else
            e_DatMatSet.hba = 0;
         end
         [ct,sigma_new,eps_new,hvar_new,aux_var] = RMapDanoSoloTraccion(eps_new,hvar_old,aux_var,...
            e_DatMatSet,e_VG);
      otherwise
         error('Modelo de da�o sol tracci�n regularizado: Modelo de ablandamiento no definido.')
   end
   
end