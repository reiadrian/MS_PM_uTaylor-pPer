% Este script permite graficar las curvas dadas en los archivos de salida
% de codigo Multiscale_Application

clear;
clc;
close all

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
% main_example_path ='/home/adrian/Documentos/Codigo_12_01/Projects.AA/Multiscale_Application/Examples/Ej22/MS_CompRVE2'; %pwd ; 
main_example_path =[pwd '/Examples/Desai/MS_RVE_Unit/'];
% main_example_path ='/home/adrian/Documentos/Codigo_12_01/Projects.AA/Multiscale_Application/Examples/MS_Comp';
% main_example_path ='/home/adrian/Documentos/Codigo_12_01/Projects.AA/Multiscale_Application/Examples/MS_Comp2';
path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
numGtitle = f_ProxString_Curva(fid);
nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nGraf = nGraf{1};

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
if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
    semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'r-');
    set(gca,'XAxisLocation','top');
else
    plot(xy(2:nSteps,1),xy(2:nSteps,2),'r-');
end

grid on
xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
titulo = strcat(ejeX,'-',ejeY);
title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
hold on
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* LOOP PARA CASOS QUE nGraf > 1                                               *
%*******************************************************************************

for iGraf = 2:nGraf
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
%Archivo "Macro" de gr�ficos a ser abierto
    actGraf = num2str(iGraf,'%03d');
    macroI = strcat('Macro',actGraf,'.dat');

    main_example_path =[pwd '/Examples/Desai/Comp_Desai_SI/'];
    path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
    file_open = fullfile(path_file,file);
    fid = fopen(file_open,'r');    
       
% N�mero steps provenientes del codigo    
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
    figure(iGraf)
    if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
        semilogx(xy(:,1),xy(:,2),'r-');
        set(gca,'XAxisLocation','top');
    else
        plot(xy(:,1),xy(:,2),'r-');
    end
    
    grid on
    xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
    ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
    titulo = strcat(ejeX,'-',ejeY);
    title(titulo,'FontSize',14,'FontWeight','bold','Color','k');

%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

end

%% Seccion aplicable para comparar con curva de datos
% %*******************************************************************************
% %* SE CIERRA ARCHIVO                                                           *
% %*******************************************************************************

% tidespor = f_CurvaDato;
% figure(1);
% hold on;
% semilogx(tidespor(:,1),tidespor(:,2),'b-');
% figure(2);
% hold on;
% semilogx(tidespor(:,1),tidespor(:,3),'b-');
%%
