clear;
clc;
close all;

pos2 = 1:51;
pos = [1:11 13:20 22:29 32:38 41:47 49:51];
%*******************************************************************************
%* PORO PRESIONES                                                             *
%*******************************************************************************
%*******************************************************************************
%* %RVE 0.1x0.1   NUEVO                                                              *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Het_Case/RVE_01/'];
path_file= main_example_path; file = 'Macro004.dat' ; %'\Examples\MS_Comp\'
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
h(1,1)=semilogx(xy(pos,1),xy(pos,2),'r--o','LineWidth',2,'MarkerSize',8);
ax = set(gca,'XAxisLocation','top','FontSize',24);
grid on
ylabel('Pore pressure [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
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
%* PORO PRESIONES                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE 0.5X0.5 - NUEVO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Het_Case/RVE_05/'];
path_file= main_example_path; file = 'Macro004.dat' ; %'\Examples\MS_Comp\'
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
h(3,1)=semilogx(xy(pos,1),xy(pos,2),'r--s','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Pore pressure [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* PORO PRESIONES                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE 1X1 - NUEVO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Het_Case/RVE_05/'];
path_file= main_example_path; file = 'Macro004.dat' ; %'\Examples\MS_Comp\'
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
h(4,1)=semilogx(xy(pos,1),xy(pos,2),'r--*','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Pore pressure [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);
% 
%*******************************************************************************
%* PORO PRESIONES                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE 0.1X0.1 - VIEJO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Het_Case/RVE_01/'];
path_file= main_example_path; file = 'Macro004.dat' ; %'\Examples\MS_Comp\'
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
h(7,1)=semilogx(xy(pos,1),xy(pos,2),'b--o','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Pore pressure [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);
% 
%*******************************************************************************
%* PORO PRESIONES                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE 0.5X0.5 - VIEJO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Het_Case/RVE_05/'];
path_file= main_example_path; file = 'Macro004.dat' ; %'\Examples\MS_Comp\'
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
h(8,1)=semilogx(xy(pos,1),xy(pos,2),'b--s','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Pore pressure [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* PORO PRESIONES                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE 1X1 - VIEJO                                                         *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Het_Case/RVE_1/'];
path_file= main_example_path; file = 'Macro004.dat' ; %'\Examples\MS_Comp\'
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
h(9,1)=semilogx(xy(pos,1),xy(pos,2),'b--*','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Pore pressure [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
lgd1 = legend({'RVE (0.1x0.1)-COE','RVE (0.5x0.5)-COE','RVE (1x1)-COE',...
    'RVE (0.1x0.1)-FOE','RVE (0.5x0.5)-FOE','RVE (1X1)-FOE'},...
'Location','southeast');
lgd1.FontSize = 24;
fclose(fid);
% 
% %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% 
%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                               *
%*******************************************************************************
%*******************************************************************************
%* %RVE 0.1x0.1  NUEVO                                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Het_Case/RVE_01/'];
path_file= main_example_path; file = 'Macro032.dat' ; %'\Examples\MS_Comp\'
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
figure(2);
h2= gobjects(10,1); 
% subplot(2,1,1)
h2(1,1)=semilogx(xy(pos,1),xy(pos,2),'r--o','LineWidth',2,'MarkerSize',8);
ax = set(gca,'XAxisLocation','top','FontSize',24);
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
% ax.FontSize = 24;
% ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);


%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                               *
%*******************************************************************************
%*******************************************************************************
%* RVE 0.5X0.5 - NUEVO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Het_Case/RVE_05/'];
path_file= main_example_path; file = 'Macro032.dat' ; %'\Examples\MS_Comp\'
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
h2(3,1)=semilogx(xy(pos,1),xy(pos,2),'r--s','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                               *
%*******************************************************************************
%*******************************************************************************
%* RVE 1X1 - NUEVO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Het_CaseRVE_05/'];
path_file= main_example_path; file = 'Macro032.dat' ; %'\Examples\MS_Comp\'
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
h2(4,1)=semilogx(xy(pos,1),xy(pos,2),'r--*','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);
% 
%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                             *
%*******************************************************************************
%*******************************************************************************
%* RVE 0.1X0.1 - VIEJO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Het_Case/RVE_01/'];
path_file= main_example_path; file = 'Macro032.dat' ; %'\Examples\MS_Comp\'
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
h2(7,1)=semilogx(xy(pos,1),xy(pos,2),'b--*','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);
% 
% 
%*******************************************************************************
%* VELOCIDAD DE FILTRACION                                                              *
%*******************************************************************************
%*******************************************************************************
%* RVE_0.5X0.5 - VIEJO                                                             *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Het_Case/RVE_05/'];
path_file= main_example_path; file = 'Macro032.dat' ; %'\Examples\MS_Comp\'
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
h2(8,1)=semilogx(xy(pos,1),xy(pos,2),'b--s','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
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
%* RVE 1X1 - VIEJO                                                         *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Het_Case/RVE_1/'];
path_file= main_example_path; file = 'Macro032.dat' ; %'\Examples\MS_Comp\'
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
h2(9,1)=semilogx(xy(pos,1),xy(pos,2),'b--o','LineWidth',2,'MarkerSize',8);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time ','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
lgd1 = legend({'RVE (0.1x0.1)-COE','RVE (0.5x0.5)-COE','RVE (1x1)-COE',...
    'RVE (0.1x0.1)-FOE','RVE (0.5x0.5)-FOE','RVE (1X1)-FOE'},...
'Location','southeast');
lgd1.FontSize = 24;
fclose(fid);