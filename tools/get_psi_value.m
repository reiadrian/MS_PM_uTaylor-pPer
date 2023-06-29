function psi_value = get_psi_value(funbc,time,e_VG)

   if e_VG.protype==0 %AA
   %Dtime = e_VG.Dtime;
   %max_time = e_VG.max_time;

   %Se esta considerando que el time ingresado siempre es mayor que cero, por lo que no se hace
   %ninguna verificacion si no se cumple con esto. Tambien se asume que esta ordenados los tiempos
   %ingresados, pero si no siempre se usa el primer numero mayor (segun estan ordenados).
   npfun = size(funbc,1);

   % Cuando time es menor a la primera abscisa (tiempo) de la funcion.
   if time<funbc(1,1)
      %Puede devolver valores de Psi menores que cero (ver que hacer con esto).
      warning(['Analisis no lineal: Determinacion del factor Psi: El tiempo de calculo es menor',...
         ' al primer punto de la funcion Psi definida: Se utiliza una extrapolacion de',...
         ' la primer recta ingresada.']) %#ok<WNTAG>
   end
   % Cuando time esta dentro de los valores de la funcion definida.
   for i = 2:npfun
      %Esta forma realizar permite que se pueda ingresar funciones con saltos en Psi, solo se debe
      %ingresar para el mismo tiempo dos valores y que los mismos ordenados de menor a mayor
      %(siempre se busca el primer tiempo mayor a time, y se interpola entre este y el anterior, y
      %luego se sale de la funcion).
      if time<=funbc(i,1)
         %Interpolacion lineal entre el valor anterior (i-1) y el siguiente.
         psi_value = funbc(i-1,2)+(time-(funbc(i-1,1)))*(funbc(i,2)-funbc(i-1,2))/...
            (funbc(i,1)-funbc(i-1,1));
         %Cuando se encuentra el primer tiempo valor mayor o igual a time, y se sale (esto permite 
         %imponer funciones con salto y que no se sobreescriba el valor de Psi).
         return
      end
   end
   % Cuando time es mayor a los valores de la funcion definida.
   %Por el return anterior aca se pasa solo si ocurre que el time es mayor.
   %if time>funbc(npfun,1)
   psi_value = funbc(npfun-1,2)+(time-(funbc(npfun-1,1)))*(funbc(npfun,2)-funbc(npfun-1,2))/...
      (funbc(npfun,1)-funbc(npfun-1,1));
   warning(['Analisis no lineal: Determinacion del factor Psi: Se supera el tiempo hasta ',...
      'donde se definio la funcion: Se utiliza una extrapolacion de la ultima recta ingresada.']) %#ok<WNTAG>
   %end
   %#########################################################################
   elseif e_VG.protype==1 %AA
         if e_VG.istep==0
             psi_value = 0.0;
         else
             psi_value = funbc(e_VG.istep,1); 
         end    %AA
 %#########################################################################
   end
