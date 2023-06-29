function tipoDat = f_tipoDat(seccion,fId)
% PUEDO USAR findstr(a,'_')

switch seccion
   case 'Tiempo'                     
      tipoDat = seccion;
   case 'Componente_X_del_desplazamiento_del_nodo'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Componente X del desplazamiento del nodo ',num2str(tipoDat1)];
   case 'Componente_Y_del_desplazamiento_del_nodo'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Componente Y del desplazamiento del nodo ',num2str(tipoDat1)];
   case 'Componente_X_de_la_fuerza_interna_del_nodo'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Componente X de la fuerza del nodo ',num2str(tipoDat1)];
   case 'Componente_Y_de_la_fuerza_interna_del_nodo'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Componente Y de la fuerza del nodo ',num2str(tipoDat1)];
   case 'Componente_X_del_salto_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Componente X del salto del elemento ',num2str(tipoDat1)];
   case 'Componente_Y_del_salto_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Componente Y del salto del elemento ',num2str(tipoDat1)];
   case 'Componente_X_de_la_tracci�n_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Componente X de la tracci�n del elemento ',num2str(tipoDat1)];
   case 'Componente_Y_de_la_tracci�n_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Componente Y de la tracci�n del elemento ',num2str(tipoDat1)];
   case 'Deformacion_Exx_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Deformacion Exx del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Deformacion_Eyy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Deformacion Eyy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Deformacion_Ezz_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Deformacion Ezz del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Deformacion_Exy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Deformacion Exy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Deformacion_Eyx_del_elemento'     %Para el caso de tensores no sim�tricos (LD por ejemplo)
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId)); 
      tipoDat = ['Deformacion Eyx del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Deformacion_fluctuante_Exx_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Deformacion fluctuante Exx del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Deformacion_fluctuante_Eyy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Deformacion fluctuante Eyy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Deformacion_fluctuante_Ezz_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Deformacion fluctuante Ezz del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Deformacion_fluctuante_Exy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Deformacion fluctuante Exy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_Txx_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension Txx del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_Tyy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension Tyy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_Tzz_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension Tzz del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_Txy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension Txy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_Tyx_del_elemento'     %Para el caso de tensores no sim�tricos (LD por ejemplo)
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension Tyx del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
    % AA Tensiones efectivas
    case 'Tension_efectiva_TExx_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension efectiva TExx del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_efectiva_TEyy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension efectiva TEyy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_efectiva_TEzz_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension efectiva TEzz del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_efectiva_TExy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension efectiva TExy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
    % AA
    % AA Tensiones totales
    case 'Tension_total_TTxx_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension total TExx del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_total_TTyy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension total TEyy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_total_TTzz_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension total TEzz del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Tension_total_TTxy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Tension total TExy del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
    % AA
    % AA Poropresiones
    case 'Propresion_del_nodo'
      tipoDat1 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Propresion del nodo ',num2str(tipoDat1)];
   % AA
   % AA Logaritmo del tiempo
   case 'logT'                     
      tipoDat = 'logT';
   case 'VelocidadFiltracion_Vx_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Velocidad de Filtracion Vx del_elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'VelocidadFiltracion_Vy_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Velocidad de Filtracion Vy del_elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
      case 'VelocidadFiltracion_Vy_n+theta_Vx_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Velocidad de Filtracion Vx_n+theta del_elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'VelocidadFiltracion_Vy_n+theta_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Velocidad de Filtracion Vy_n+theta del_elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Masa_fluido_chi_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Chi ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Variable_de_da�o_del_elemento'
      tipoDat1 = str2double(f_ProxString_Curva(fId)); 
      f_ProxString_Curva(fId);
      tipoDat3 = str2double(f_ProxString_Curva(fId));
      tipoDat = ['Variable de da�o del elemento ',num2str(tipoDat1),' y del punto de gauss ','',num2str(tipoDat3)];
   case 'Disipacion_estructural_del_test'
      tipoDat = 'Disipacion estructural del test';
    case 'Homogenized_Strain_(component_xx)'
      tipoDat = 'Homogenized Strain (component xx)';
    case 'Homogenized_Strain_(component_yy)'
      tipoDat = 'Homogenized Strain (component yy)';
    case 'Homogenized_Strain_(component_zz)'
      tipoDat = 'Homogenized Strain (component zz)';
    case 'Homogenized_Strain_(component_xy)'
      tipoDat = 'Homogenized Strain (component xy)';
    case 'Homogenized_Stress_(component_xx)'
      tipoDat = 'Homogenized Stress (componentxx)';
    case 'Homogenized_Stress_(component_yy)'
      tipoDat = 'Homogenized Stress (component yy)';
    case 'Homogenized_Stress_(component_zz)'
      tipoDat = 'Homogenized Stress (component zz)';
    case 'Homogenized_Stress_(component_xy)'
      tipoDat = 'Homogenized Stress (component xy)';
    case 'Homogenized_Stress_(component_nn)'
      tipoDat = 'Homogenized Stress (component nn)';
    case 'Homogenized_Stress_(component_nt)'
      tipoDat = 'Homogenized Stress (component nt)';
    case 'Mean_spherical_stress'
      tipoDat = 'Mean spherical stress';
    case 'J2'
      tipoDat = 'J2';
   otherwise
      error('Archivos de datos: Inicializaci�n: No est� definido este tipo de dato.')
end