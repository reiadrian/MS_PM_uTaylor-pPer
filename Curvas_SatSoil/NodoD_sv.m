clear;
clc;
close all;

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Monoescala/'];
path_file= main_example_path; file = 'Macro027.dat' ; %'\Examples\MS_Comp\'
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
h= gobjects(5,2); 
% subplot(2,1,1)
h(1,1)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'b-','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Effective stress 22 [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Rel_Const_Actual/RVE1x1_UCR/'];
path_file= main_example_path; file = 'Macro027.dat' ; %'\Examples\MS_Comp\'
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
h(2,1)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'r-','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Effective stress 22 [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Rel_Const_Actual/RVE01x01_UCR/'];
path_file= main_example_path; file = 'Macro027.dat' ; %'\Examples\MS_Comp\'
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
h(3,1)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'r--','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Effective stress 22 [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Rel_Const_Anterior/RVE1x1_OCR/'];
path_file= main_example_path; file = 'Macro027.dat' ; %'\Examples\MS_Comp\'
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
h(4,1)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'k-','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Effective stress 22 [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Rel_Const_Anterior/RVE01x01_OCR/'];
path_file= main_example_path; file = 'Macro027.dat' ; %'\Examples\MS_Comp\'
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
h(5,1)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'k--','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Effective stress 22 [kPa]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
lgd1 = legend({'Mono-scale','RVE (1x1)-URC','RVE (0.1X0.1)-URC','RVE (1x1)-ORC','RVE (0.1X0.1)-ORC'},...
'Location','northeast');
lgd1.FontSize = 24;
fclose(fid);

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Monoescala/'];
path_file= main_example_path; file = 'Macro041.dat' ; %'\Examples\MS_Comp\'
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
% h= gobjects(5,2); 
% subplot(2,1,2)
h(1,2)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'b-','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Rel_Const_Actual/RVE1x1_UCR/'];
path_file= main_example_path; file = 'Macro041.dat' ; %'\Examples\MS_Comp\'
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
h(2,2)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'r-','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Rel_Const_Actual/RVE01x01_UCR/'];
path_file= main_example_path; file = 'Macro041.dat' ; %'\Examples\MS_Comp\'
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
h(3,2)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'r--','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Rel_Const_Anterior/RVE1x1_OCR/'];
path_file= main_example_path; file = 'Macro041.dat' ; %'\Examples\MS_Comp\'
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
h(4,2)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'k-','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Rel_Const_Anterior/RVE01x01_OCR/'];
path_file= main_example_path; file = 'Macro041.dat' ; %'\Examples\MS_Comp\'
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
h(5,2)=semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'k--','LineWidth',2);
set(gca,'XAxisLocation','top');
grid on
ylabel('Flux velocity [m/s]','FontSize',24,'FontWeight','bold','Color','k');
xlabel('Time [s]','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
lgd2 = legend({'Mono-scale','RVE (1x1)-URC','RVE (0.1X0.1)-URC','RVE (1x1)-ORC','RVE (0.1X0.1)-ORC'},...
'Location','northeast');
lgd2.FontSize = 24;
fclose(fid);
