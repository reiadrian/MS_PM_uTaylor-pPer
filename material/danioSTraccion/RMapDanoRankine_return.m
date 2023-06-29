function [ct,sigma_new,eps_new,hvar_new,aux_var] = RMapDanoRankine_return(eps_new,hvar_old,aux_var,...
   e_DatMatSet,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DE FUERZA INTERNA Y MATRIZ TANGENTE: PLANE STRAIN - 3D   *
%*  MODELO DE DA�O ESCALAR CON ABLANDAMIENTO                                              *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% Variables globales
factor = 1e-6;

% Recupera propiedades del material
E = e_DatMatSet.young;
%sigma_u = e_DatMatSet.ftult;
H = e_DatMatSet.hba;
tita = e_DatMatSet.tit;
ce = e_DatMatSet.ce;
r_0 = e_DatMatSet.r_0;
esImplex = e_DatMatSet.esImplex;

% Recupera variables internas
r_old = hvar_old(1);
q_old = hvar_old(2);

% Inicializacion de variables
%r_0 = sigma_u/sqrt(E);
q_0 = r_0;
zero_q = factor*r_0;
%Inicializaci�n de variables old.
if r_old<=0 
   r_old = r_0;
   q_old = r_0;
end

%Si se fuerza que se comporte el�sticamente, la definici�n puede venir globalmente o por puntos de gauss.
if esImplex
   elast = aux_var(4);
else
   elast = aux_var(3);
end
elast = e_VG.elast||elast;

% Tension efectiva
effective_stress = ce*eps_new;

% Tensi�n efectiva positiva
if e_VG.ndime==2
   %Para deformaci�n plana y tensi�n plana
   [m_max_Stress, tensorDirPrinMax] = f_AutoValMax2D(effective_stress,e_VG);
%    m_SigmaMatriz = [effective_stress(1),effective_stress(4);effective_stress(4),effective_stress(2)];
%    [v,d] = eig(m_SigmaMatriz);
%    d(d<0) = 0;
%    m_SigmaMatriz = v*d*v';
%    m_StressEfecPos = zeros(e_VG.ntens,1);
%    m_StressEfecPos(1) = m_SigmaMatriz(1);
%    m_StressEfecPos(2) = m_SigmaMatriz(4);
%    m_StressEfecPos(4) = m_SigmaMatriz(2);
%    if effective_stress(3)>0
%       m_StressEfecPos(3) = effective_stress(3);
%    end
else
   error('Da�o solo tracci�n: No definido los autovalores positivo en el caso 3D.')
end

% Estado Trial
inv_sqrt_E = 1/sqrt(E);
r_trial =inv_sqrt_E * m_max_Stress;

if r_trial>r_old&&~elast          
  
   % Da�o
   fload_new = 1.0; 
   
   r_new = r_trial;
   
   if (tita == 0)
      % Caso lineal
      %delta_r = r_trial-r_old; 
      %q_new = q_old+H*delta_r;
      %Para el caso lineal escribirlo en forma incremental o como la evaluaci�n de la funci�n q(r) es lo
      %mismo, y para que est� consistente como se escribi� en forma exponencial, se escribe como q(r).
      q_new = q_0+H*(r_new-r_0);
      %En el caso que en el ablandamiento lineal se llegue por debajo de zero_q se fuerza que comporte como
      %exponencial, as� sigue siendo un modelo de da�o pero con derivadas continuas y sin llegar a zero las
      %tensiones y los tensores tangentes constitutivos. 
      %Se verifica que H no sea cero (aunque no deber�a serlo en da�o). En este caso no se hace ning�n cambio.
      %La energ�a de fractura no ser�a exactamente la esperada, ya nunca llegar�a a cero la gr�fica
      %lineal.
      %Verificar que funciona bien!!
%       if q_new<zero_q&&H
%          %fprintf('El elemento %d y PG %d lleg� al l�mite de q (modelo de da�o).\n',e_VG.iElemNum,e_VG.iPG)
%          %Se determina r_0 justo para la posici�n de q_0, que servir� como punto inicial de la curva
%          %exponencial.
%          r_0 = (zero_q-q_0)/H+r_0;
%          q_0 = zero_q;
%          %Esto no funciona cuando exp(H*(r_new-r_0)/q_0) es aproximadamente exp(-745), ya que un poco mayor da
%          %cero el resultado del exponencial por redondeo num�rico (probar exp(-746), devolviendo que qnew es cero.
%          %Hay que determinar bien que factor usar para que no se pase esto.
%          q_new = q_0*exp(H*(r_new-r_0)/q_0);
%          tita = 1;
%          %Se cambia el fload para indicar esta situaci�n.
%          fload_new = 2.0;
%       end  
   elseif (tita == 1)
      % Caso exponencial
      q_new = q_0*exp(H*(r_new-r_0)/q_0);
   end
   
   %beta_1 = 1.0;
   if q_new<zero_q
      %Si el q_new baja de cierto valor se obliga al material que se comporte el�sticamente (no da�a) y tiene
      %un tensor constitutivo tangente el�stico.
      q_new = zero_q;
      %Esto determinaciones de r para zero_q falla si H=0, pero en ese caso q_new no evoluciona y nunca llega
      %a zero_q (es un factor de q_0 que es el valor inicial de q).
      if tita==0
         r_new = (zero_q-q_0)/H+r_0;
      elseif tita==1
         r_new = q_0*log(zero_q/q_0)/H+r_0;
      end
      beta_2 = q_new/r_new;
      beta_3 = 0;
      %Se cambia el fload para indicar esta situaci�n.
      fload_new = 2.0;
   else
      beta_2 = q_new/r_new; 
      if (tita == 0)
         beta_3  = (H*r_new-q_new)/r_new^2;
      elseif (tita == 1)
         beta_3  = q_new/r_new^2 * (H*r_new/q_0 - 1);
      end
   end
   
   %damage_new = 1.0 - (q_new/r_new);
   damage_new = 1.0-beta_2;
   sigma_new  = beta_2*effective_stress;
   % Modulo elastoplastico consistente
   %ct = beta_2*ce + beta_1*beta_3*(m_StressEfecPos*effective_stress');
   m_StressEfecPos = inv_sqrt_E * ce*tensorDirPrinMax;
   ct = beta_2*ce+beta_3*(effective_stress*m_StressEfecPos');
    
else
    
   % Descarga elastica
   if (r_old > r_0)
      % Esto es para descarga el�stica
      fload_new  = -1.0; 
   else
      fload_new  = 0.0; 
   end
   r_new = r_old;
   q_new = q_old;
   if E~=0
      beta_2 = q_new/r_new;
   else
      %Ac� se est� considerando que si E=0 nunca da�a. Se puede ver que r_trial>r_old nunca va a verificar ya
      %que el primer valor de r_old es r_0 = sigma_u/sqrt(E) = Inf, y rtrial = sqrt(eps_new'*0) = 0, y por lo
      %tanto nunca entra en la zona de da�o.
      beta_2 = 1;
   end
   damage_new = 1-beta_2;
   
   %Tensi�n
   sigma_new = beta_2*effective_stress;
   %beta_1     = 0.0;
   %beta_3     = 0.0;
   
   % Modulo elastoplastico consistente
   ct = beta_2*ce;

end

%Se realiza c�lculo de las variables Impl�citas-Expl�citas (Implex) si se seleccion� el c�lculo de esta forma.
if esImplex
   facDelta_r = aux_var(3);
   if ~elast
      delta_r_old = hvar_old(3);
      r_newImplex = r_old+delta_r_old*facDelta_r;
      if tita==0
         q_newImplex = q_0+H*(r_newImplex-r_0);
      elseif tita==1
         q_newImplex = q_0*exp(H*(r_newImplex-r_0)/q_0);
      end
      sigma_newImplex = q_newImplex/r_newImplex*effective_stress;
      % Se reemplaza ya que no se utiliza las sigma impl�cita por las implex en el argumento de salida.
      sigma_new = sigma_newImplex;
      ctImplex = q_newImplex/r_newImplex*ce;
      %Ver c�mo definir para el implex el damage_new y el fload_new.
      %Se almacena los resultados.
      ct = struct('Implex',ctImplex,'Impli',ct);
   else
      ct = struct('Implex',ct,'Impli',ct);
   end
   %Valores impl�citos del incremento de r
   delta_r = r_new-r_old;
   % Variables hist�ricas
   %(se guarda los valores impl�citos)
   hvar_new = [r_new;q_new;delta_r;hvar_old(4:6)];
   % Variables auxiliares
   aux_var = [damage_new;fload_new;facDelta_r;elast];
else
   %Variables hist�ricas
   hvar_new = [r_new;q_new;hvar_old(3:5)];
   % Variables auxiliares
   aux_var = [damage_new;fload_new;elast];
end

% Variables historicas
%hvar_new = [r_new ; q_new ; damage_new ; fload_new];
% Variables auxiliares
%aux_var = [damage_new;fload_new];

