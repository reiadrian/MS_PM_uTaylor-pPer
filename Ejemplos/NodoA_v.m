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
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/Pruebas/RVE_1_Taylor/'];
path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

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
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/Pruebas/RVE_1_Lineal/'];
path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

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

%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE_1X1 - NUEVO                                                         *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/Pruebas/RVE_1_Per_Tay_16/'];
path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

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
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Beta1_Homog/Pruebas/RVE_1_Periodico/'];
path_file= main_example_path; file = 'Macro034.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

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
h(7,1)=semilogx(xy(1:nSteps,1),xy(1:nSteps,2),'b--s','LineWidth',2,'MarkerSize',8);
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


lgd1 = legend({'Mono-scale','Taylor','Lineal','u-Per / p-Taylor (Numer)','Periodico'},...
'Location','northeast');
lgd1.FontSize = 24;
