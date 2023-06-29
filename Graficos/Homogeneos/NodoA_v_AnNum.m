clc;
close all;

%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                             *
%*******************************************************************************
%*******************************************************************************
%* %MONOESCALA                                                                 *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Monoescala/'];
path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
% numGtitle = f_ProxString_Curva(fid);
% nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% nGraf = nGraf{1};

% Numero steps provenientes del codigo
numSteps = f_ProxString_Curva(fid);
nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nSteps = nSteps{1}; 

% Lectura del t�tulo para eje de abscisas
seccion = f_ProxString_Curva(fid);
ejeX = f_tipoDat(seccion,fid);
vs = f_ProxString_Curva(fid);

% Lectura del t�tulo para eje de aordenadas
seccion = f_ProxString_Curva(fid);
ejeY = f_tipoDat(seccion,fid);

% Lectura de datos XY
format = '%q %q';
xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
xy = xy{1};
xy = str2double(xy(:,:));
% Grafico 1
figure(1);
h= gobjects(10,1); 
% subplot(2,1,1)
h(1,1)=semilogx(xy(1:nSteps,1),xy(1:nSteps,2),'k-x','LineWidth',2,'MarkerSize',8);
ax = set(gca,'XAxisLocation','top','FontSize',24);
grid on
ylabel('Flux velocity [m/day]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [day] ','FontSize',24,'FontWeight','bold','Color','k');
% ax.FontSize = 24;
% ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%###############################################################################
%###############################################################################

%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE_0.1X0.1 - NUEVO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/An_Num/Analitico/RVE_10_Analitico/'];
path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
% numGtitle = f_ProxString_Curva(fid);
% nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% nGraf = nGraf{1};

% Numero steps provenientes del codigo
numSteps = f_ProxString_Curva(fid);
nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nSteps = nSteps{1}; 

% Lectura del t�tulo para eje de abscisas
seccion = f_ProxString_Curva(fid);
ejeX = f_tipoDat(seccion,fid);
vs = f_ProxString_Curva(fid);

% Lectura del t�tulo para eje de aordenadas
seccion = f_ProxString_Curva(fid);
ejeY = f_tipoDat(seccion,fid);

% Lectura de datos XY
format = '%q %q';
xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
xy = xy{1};
xy = str2double(xy(:,:));
% Grafico 1
h(3,1)=semilogx(xy(1:nSteps,1),xy(1:nSteps,2),'r--*','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/day]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [day] ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);



%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE_05X05 - NUEVO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/An_Num/Analitico/RVE_10_Analitico/'];
path_file= main_example_path; file = 'Macro050.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
% numGtitle = f_ProxString_Curva(fid);
% nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% nGraf = nGraf{1};

% Numero steps provenientes del codigo
numSteps = f_ProxString_Curva(fid);
nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nSteps = nSteps{1}; 

% Lectura del t�tulo para eje de abscisas
seccion = f_ProxString_Curva(fid);
ejeX = f_tipoDat(seccion,fid);
vs = f_ProxString_Curva(fid);

% Lectura del t�tulo para eje de aordenadas
seccion = f_ProxString_Curva(fid);
ejeY = f_tipoDat(seccion,fid);

% Lectura de datos XY
format = '%q %q';
xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
xy = xy{1};
xy = str2double(xy(:,:));
% Grafico 1
h(4,1)=semilogx(xy(1:nSteps,1),xy(1:nSteps,2),'b--*','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/day]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [day] ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);
% lgd1 = legend({'V^{n+1} - Mono-scale','V^{n+1} - Analitico','V^{n+\theta} - Analitico'},...
% 'Location','northeast');
% lgd1.FontSize = 24;
% lgd1 = legend({'V^{n+1} - Mono-scale','V^{n+1} - Analitico','V^{n+1} - Numerico'},...
% 'Location','northeast');
% lgd1.FontSize = 24;


%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE_1X1 - NUEVO                                                         *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/An_Num/Numerico/RVE_10_Numerico/'];
path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
% numGtitle = f_ProxString_Curva(fid);
% nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% nGraf = nGraf{1};

% Numero steps provenientes del codigo
numSteps = f_ProxString_Curva(fid);
nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nSteps = nSteps{1}; 

% Lectura del t�tulo para eje de abscisas
seccion = f_ProxString_Curva(fid);
ejeX = f_tipoDat(seccion,fid);
vs = f_ProxString_Curva(fid);

% Lectura del t�tulo para eje de aordenadas
seccion = f_ProxString_Curva(fid);
ejeY = f_tipoDat(seccion,fid);

% Lectura de datos XY
format = '%q %q';
xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
xy = xy{1};
xy = str2double(xy(:,:));
% Grafico 1
h(5,1)=semilogx(xy(1:nSteps,1),xy(1:nSteps,2),'r--o','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/day]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [day] ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%###############################################################################
%###############################################################################

%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE_0.1X0.1 - VIEJO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/An_Num/Numerico/RVE_10_Numerico/'];
path_file= main_example_path; file = 'Macro050.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
% numGtitle = f_ProxString_Curva(fid);
% nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% nGraf = nGraf{1};

% Numero steps provenientes del codigo
numSteps = f_ProxString_Curva(fid);
nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nSteps = nSteps{1}; 

% Lectura del t�tulo para eje de abscisas
seccion = f_ProxString_Curva(fid);
ejeX = f_tipoDat(seccion,fid);
vs = f_ProxString_Curva(fid);

% Lectura del t�tulo para eje de aordenadas
seccion = f_ProxString_Curva(fid);
ejeY = f_tipoDat(seccion,fid);

% Lectura de datos XY
format = '%q %q';
xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
xy = xy{1};
xy = str2double(xy(:,:));
% Grafico 1
h(7,1)=semilogx(xy(1:nSteps,1),xy(1:nSteps,2),'b--o','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/day]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [day] ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);
% 
lgd1 = legend({'Mono-scale','V^{n+1} - Analitico','V^{n+\theta} - Analitico','V^{n+1} - Numerico','V^{n+\theta} - Numerico'},...
'Location','northeast');
lgd1.FontSize = 24;
% lgd1 = legend({'Mono-scale','RVE 1x1-V^{n+\theta}','RVE 1x1-V^{n+1}','RVE 0.5x1-V^{n+\theta}','RVE 0.5x0.5-V^{n+1}'},...
% 'Location','northeast');
% lgd1.FontSize = 24;
% lgd1 = legend({'Mono-scale','V^{n+\theta}=f(Chi*y/Dtime)','V^{n+1}=f(Chi*y)','V^{n+1}=f(Chi*y/Dtime)','V^{n+1}=f(V^{n+theta})'},...
% 'Location','northeast');
% lgd1.FontSize = 24;
% 
% 
% 
% 
% 
% 
% 

% 
%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE_0.5X0.5 - VIEJO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/ntheta/Monoescala/'];
path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
% numGtitle = f_ProxString_Curva(fid);
% nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% nGraf = nGraf{1};

% Numero steps provenientes del codigo
numSteps = f_ProxString_Curva(fid);
nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nSteps = nSteps{1}; 

% Lectura del t�tulo para eje de abscisas
seccion = f_ProxString_Curva(fid);
ejeX = f_tipoDat(seccion,fid);
vs = f_ProxString_Curva(fid);

% Lectura del t�tulo para eje de aordenadas
seccion = f_ProxString_Curva(fid);
ejeY = f_tipoDat(seccion,fid);

% Lectura de datos XY
format = '%q %q';
xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
xy = xy{1};
xy = str2double(xy(:,:));
% Grafico 1
h(8,1)=semilogx(xy(1:nSteps,1),xy(1:nSteps,2),'b--s','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/day]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [day] ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);
% % 
% % %*******************************************************************************
% % %* VELOCIDAD DE FILTRACION                                                             *
% % %*******************************************************************************
% % %*******************************************************************************
% % %* RVE_1X1 - VIEJO                                                         *
% % %*******************************************************************************
% % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % DIRECTORIO EN LA LINEA INFERIOR
% % main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/RVE_05_Homog_FOE (otra copia)/'];
% % path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
% % file_open = fullfile(path_file,file);
% % fid = fopen(file_open,'r');
% % 
% % % Numero de graficos provenientes del codigo
% % % numGtitle = f_ProxString_Curva(fid);
% % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nGraf = nGraf{1};
% % 
% % % Numero steps provenientes del codigo
% % numSteps = f_ProxString_Curva(fid);
% % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % nSteps = nSteps{1}; 
% % 
% % % Lectura del t�tulo para eje de abscisas
% % seccion = f_ProxString_Curva(fid);
% % ejeX = f_tipoDat(seccion,fid);
% % vs = f_ProxString_Curva(fid);
% % 
% % % Lectura del t�tulo para eje de aordenadas
% % seccion = f_ProxString_Curva(fid);
% % ejeY = f_tipoDat(seccion,fid);
% % 
% % % Lectura de datos XY
% % format = '%q %q';
% % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % xy = xy{1};
% % xy = str2double(xy(:,:));
% % % Grafico 1
% % h(9,1)=semilogx(xy(1:nSteps,1),xy(1:nSteps,2),'b--o','LineWidth',2,'MarkerSize',8);
% % set(gca,'XAxisLocation','top');
% % grid on
% % ylabel('Flux velocity [m/day]','FontSize',24,'FontWeight','bold','Color','k');
% % xlabel('Time [day] ','FontSize',24,'FontWeight','bold','Color','k');
% % ax.FontSize = 24;
% % ay.FontSize = 24;
% % hold on
% % %*******************************************************************************
% % %* SE CIERRA ARCHIVO                                                           *
% % %*******************************************************************************
% % lgd1 = legend({'Mono-scale','RVE (0.1x0.1)-SOE','RVE (0.5x0.5)-SOE',...
% %     'RVE (1X1)-SOE','RVE (0.1x0.1)-FOE','RVE (0.5x0.5)-FOE','RVE (1X1)-FOE'},...
% % 'Location','northeast');
% % lgd1.FontSize = 24;
% % fclose(fid);
% % % 
% % % % %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% % % % %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% % % % %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&