function [PGLOC,PGNLOC,PGLOCT] = get_pgloc_pgnloc(hvar_old,e_VG)

   npg   = e_VG.npg;
   ntens = e_VG.ntens;

   PGLOC            = get_fload(hvar_old,e_VG);
   
   %Se asume que la funci�n get_fload devuelve los PG localizados, por lo que se fuerza a que si un
   %PG del elemento est� localizado, todos los PG del mismo est�n como localizado (es decir elemento
   %localizado).
   PGLOC = repmat(any(PGLOC,1),npg,1);
   PGNLOC = ~PGLOC;
   %Considerando las cantidad de componentes de tensi�n y deformaci�n.
   PGLOCT = PGLOC(repmat(1:npg,ntens,1),:);
   
end