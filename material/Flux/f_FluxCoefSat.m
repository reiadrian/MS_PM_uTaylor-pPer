function e_DatSet = f_FluxCoefSat(e_DatSet)
% Calcula el valor de la constatne alpha de Biot y la constante que toma en
% cuenta la compresibilidad del solido y del fluido

 e_DatMatSet = e_DatSet.e_DatMat;
 e0 = e_DatMatSet.e0; % Relacon de vacions inicial
 Ks = e_DatMatSet.Ks; % Compresibilidad del solido(grano)/Rigidez del grano
 Kw = e_DatMatSet.Kw; % Compresibilidad del agua
 kx = e_DatMatSet.kx; % Permeabilidad en la direccion "x"
 ky = e_DatMatSet.ky; % Permeabilidad en la direccion "y"
 ce = e_DatMatSet.ce; % Matriz constitutiva elastica
 
 n = e0/(e0+1); % Porosidad
 Satgr = 1.0; %Grado de saturacion (Adopto igual a 1...solo resuelvo casos saturados)
 m = [1 1 1 0]; % Vectoriza para tensiones xx,yy,zz,xy

 DRCom = m*ce*m' ;
 DRCom = DRCom/9 ;
 Lvalu = 1.0 - DRCom/Ks ;
 
 alpha = Satgr*Lvalu ;
 
 e_DatSet.e_DatMat.m_Biot = m'.*alpha ; % Matriz de Biot
 e_DatSet.e_DatMat.beta = (Lvalu-n)*Satgr/Ks + n*Satgr/Kw - DRCom/(3*Ks^2) ; % "Coeficiente de compresibilidad". Adicione el ultimo termino conforme TESIS Di Rado
 e_DatSet.e_DatMat.m_PermK = [kx 0.0; 0.0 ky]; % Matriz de permeabilidad
 