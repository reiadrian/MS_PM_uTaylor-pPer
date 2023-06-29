function [istep_bif,m_thetaBif,m_mindetQ] = get_bifurcation_condition(m_CTHomog,istep_bif,e_VG)

delTheta = 0.0025*pi;

%ntens = e_VG.ntens;
%ndime = e_VG.ndime;
theta = (-pi/2):delTheta:(pi/2-delTheta/2);
nTheta = size(theta,2);
detQ = zeros(1,nTheta);
% n_m   = zeros(ntens,ndime);

%Versi�n semi-vectorizada.
m_nVoigt = f_Vec2Voigt2DVect([cos(theta)',sin(theta)'],e_VG);
for ind = 1:nTheta
   %Tensor de localizacion
   detQ(ind) = det(m_nVoigt(:,:,ind)'*m_CTHomog*m_nVoigt(:,:,ind));
end

% hold on
% plot(theta,detQ)
% grid on
nomArch = [e_VG.fileCompleto,'.pasBif'];
fId = fopen(nomArch,'at');
fprintf(fId,'#Paso: %d\n',e_VG.istep);
fprintf(fId,'#�ngulo Det(Q)\n');
fprintf(fId,'%f %f\n',[theta;detQ]);
fprintf(fId,'\n\n');
fclose(fId);

m_thetaBif = [];
m_mindetQ = [];
%
[mindetQ1,m_IndBif1] = min(detQ);
%
if (mindetQ1 < 0 && istep_bif == 0)
   istep_bif = e_VG.istep;
   %Se utiliza un detQ auxiliar para determinar los pr�ximos m�nimos, donde los m�nimos encontrados
   %se los reeamplaza con infinito, y as� ignorarlos. Y el detQ original se utiliza para calcular
   %las pendientes (si no se hace esto, al reemplazar con inf, se puede detectar m�nimos en forma
   %err�nea).
   detQAux = detQ;
   detQAux(m_IndBif1) = inf;
   sdoMin = 0;
   %Se hace la b�squeda del segundo m�nimo. Itera hasta encontrarlo (podr�a llegar a encontrar un
   %m�nimo local, pero supuestamente este no existe). Notar que ignora los m�ximos locales por el 
   %control de las derivadas.
   for i = 1:nTheta-1
      [mindetQ2,m_IndBif2] = min(detQAux);
      %Para considerar el anterior y el siguiente del primero y el �ltimo de la matriz, se ingresa
      %condiciones adicionales
      if (m_IndBif2==1&&detQ(1)-detQ(nTheta)<0&&detQ(2)-detQ(1)>0)||...
            (m_IndBif2==nTheta&&detQ(nTheta)-detQ(nTheta-1)<0&&detQ(1)-detQ(nTheta)>0)||...
            (m_IndBif2>1&&m_IndBif2<nTheta&&detQ(m_IndBif2)-detQ(m_IndBif2-1)<0&&...
            detQ(m_IndBif2+1)-detQ(m_IndBif2)>0) 
         sdoMin = 1;
         break
      end
      detQAux(m_IndBif2) = inf;
   end
   if ~sdoMin
      mindetQ2 = [];
      m_IndBif2 = [];
   end
   m_thetaBif = theta([m_IndBif1,m_IndBif2]);
   m_mindetQ = [mindetQ1,mindetQ2];
   %Imprensi�n en pantalla
   fprintf('An�lisis de bifurcaci�n:\n')
   fprintf('theta=%f Det(Q)=%f\n',[m_thetaBif;m_mindetQ])
   %Imprensi�n en archivo de paso que bifurca
   nomArch = [e_VG.fileCompleto,'.pasBif'];
   fId = fopen(nomArch,'at');
   %Se # porque es el comentario del GNUPlot.
   fprintf(fId,'#################################################\n');
   fprintf(fId,'#Bifurcaci�n, paso: %d (t=%f)\n',e_VG.istep,e_VG.Dtime*e_VG.istep);
   fprintf(fId,'#theta=%f Det(Q)=%f\n',[m_thetaBif;m_mindetQ]);
   fprintf(fId,'#################################################\n');
   fclose(fId);
   %plot(theta,detQ)
   %pause
end

end