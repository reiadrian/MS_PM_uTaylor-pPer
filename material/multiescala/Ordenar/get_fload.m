function PGLOC = get_fload(hvar_old,e_VG)

   conshyp = e_VG.conshyp;

   switch conshyp
      case 2
         PGLOC = get_fload_plasJ2(hvar_old,e_VG);
      case {4,10}
         PGLOC = get_fload_damage(hvar_old,e_VG);
      case 8
         PGLOC = f_GetLoad_CentralForces(hvar_old,e_VG);
      case 9
         PGLOC = get_fload_plasJ2bl(hvar_old,e_VG);
      otherwise
         error(['An�lisis NoLineal: Determinaci�n de elementos localizados: ',...
            'Modelo constitutivo no definido'])
   end
   
end
