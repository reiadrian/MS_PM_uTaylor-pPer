function m_sigmaHomog = f_HomogFront(Fint,xx,area,e_VG)

   ndime = e_VG.ndime;
   struhyp = e_VG.struhyp;
   %Hasta que se pase los nodos de la frontera, se calcula los para la Fint de todos los Nodos,
   %asumiento que las fuerzas internas en nodos interiores siempre es nula.
   %Es necesario pasarla como una matriz full para redimensionarla en tercera dimensi�n.
   m_sigmaHomog = sum(bsxfun(@times,reshape(full(Fint),ndime,1,[]),permute(xx(:,1:ndime),...
      [3,2,1])),3);
   m_sigmaHomog = 0.5*(m_sigmaHomog+m_sigmaHomog')/area;
   % Se escribe en notaci�n de Voight
   switch struhyp
      case 1
         %El vector tensor en una frontera con normal en el plano x-y no tiene componente
         %en direcci�n z, al ser deformaci�n plana (hay simetr�a en direcci�n z en 
         %cualquier plano xy que se elija). Por ello su valor es nulo.
         m_sigmaHomog = [m_sigmaHomog(1);m_sigmaHomog(4);0;m_sigmaHomog(2)];
      otherwise
         error('An�lisis no lineal: Homgenizaci�n por valores de frontera: Hip�tesis o estado de c�lculo no definido.')
   end
   
end