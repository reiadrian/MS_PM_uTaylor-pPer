function  m_VecVoigt = f_Vec2Voigt2D(m_Vector,e_VG)

   % Vector en formato matricial, asumiendo notaci�n de Voigt con el que el tensor Sigma y
   % Deformaciones teniendo la siguiente forma [Sxx;Syy;Szz;Sxy] (caso 2D, de tensi�n y deformaci�n
   % plana).

   m_VecVoigt = zeros(e_VG.ntens,e_VG.ndime);
   m_VecVoigt(1,1) = m_Vector(1);  
   m_VecVoigt(4,1) = m_Vector(2);
   m_VecVoigt(2,2) = m_Vector(2); 
   m_VecVoigt(4,2) = m_Vector(1);  

end