function m_hReg = f_LongReg(hReg,nomArchHReg,dirDat,m_NumElemSet,nElemSet)

   %Se utiliza NaN para verificar que todos los elementos tiene indicados como proceder para el c�lculo de la
   %distancia de regularizaci�n.
   %autom�ticamente, seg�n alg�n criterio posterior y impuesto en el archivo de datos.
   m_hReg = nan(1,nElemSet);
   if ~isnan(hReg)
      %En hReg se puede indicar la distancia de regularizaci�n, � en valor negativo el m�todo a usar para el
      %c�lculo autom�tico. Este valor se aplica a todos los elementos.
      m_hReg(:) = hReg;
   end
   %En la propiedad ARCHHREG se puede indicar "", quedando la variable como un string vac�o.
   if ~isempty(nomArchHReg)&&~any(isnan(nomArchHReg))
      %Se lee un archivo que se indica las distancias de regularizaci�n. Esta sobreescriben las indicadas por
      %el hReg.
      %Si se ingresa el nombre del archivo a leer (tiene que tener extensi�n hReg).
      fiD = fopen(fullfile(dirDat,[nomArchHReg,'.hReg']));
      %En el archivo se debe indicar el n�mero de elemento (el n�mero de elemento debe corresponder
      %al utilizado en el archivo de datos, sin importar su posici�n en la lista de conectividades)
      %y el distancia de regularizaci�n (una fila por cada n�mero de elemento finito).
      %No es necesario crear un archivo por cada set si se tiene varios set de da�o en la misma
      %estructura, ya los elementos se descartan si no pertenece al set.
      m_ElemLeido = fscanf(fiD,'%f %f\n',[2,inf]);
      fclose(fiD);
      %Se obtiene la posici�n en la lista de elementos del set de los elementos (denominaci�n) 
      %indicados.
      %Es transformar como se denomina los elementos en el archivo de datos (una lista no
      %necesariamente ordenada y completa, pero �nica) en su �ndice posici�n dentro de la lista de 
      %elementos del set (que es lo que se utiliza dentro del programa para indentificar al
      %elemento).
      m_IndElem = bsxfun(@eq,m_NumElemSet',m_ElemLeido(1,:));
      [m_NroIndhReg,~] = find(m_IndElem);
      %Se reemplaza �nicamente los elementos que se indic� en el archivo ARCHHREG, en el resto se impone el
 
     %valor de HREG o por defecto. Cuidado no tira ninguna notificaci�n si en el archivo ARCHHREG se olvida
      %de ingresar alg�n elemento.
      m_hReg(m_NroIndhReg) = m_ElemLeido(2,any(m_IndElem,1));
   end
   %
   if any(isnan(m_hReg))
      error(['Lectura de datos: Propiedades Materiales: Modelo constitutivo ',...
         'regularizado: Longitud de regularizaci�n: Hay elementos finitos que no tiene definido la distancia.'])
   end

end