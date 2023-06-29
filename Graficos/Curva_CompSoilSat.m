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
main_example_path =[pwd '/Examples/Desai/Desai_log/Comp_Desai_SI/'];
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
% xy(:,2) = xy(:,2)/0.3048;
% Grafico 1
figure(1);
if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
    semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'r-*');
    set(gca,'XAxisLocation','top');
else
    plot(xy(2:nSteps,1),xy(2:nSteps,2),'r-*');
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
    
%     main_example_path ='/home/adrian/Documentos/Codigo_12_01/Projects.AA/Multiscale_Application/Examples/Ej22/MS_CompRVE2'; %pwd ;
    main_example_path =[pwd '/Examples/Desai/Desai_log/Comp_Desai_SI/'];
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
%     switch iGraf
%         case {3,5,7,9,11,13,15}
%             xy(:,2) = xy(:,2)/0.3048;
%         case {2,4,6,8,10,12,14,16}
%             xy(:,2) = 1000*xy(:,2)/(47.8803);
%     end
    % Grafico 1
    figure(iGraf)
    if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
        semilogx(xy(:,1),xy(:,2),'r-*');
        set(gca,'XAxisLocation','top');
    else
        plot(xy(:,1),xy(:,2),'r-*');
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
% 
% tidespor = f_CurvaDato;
% figure(10);
% hold on;
% semilogx(tidespor(:,1),tidespor(:,2),'b-');
% figure(11);
% hold on;
% semilogx(tidespor(:,1),tidespor(:,3),'b-');
%%
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =[pwd '/Examples/Desai/Desai_log/Comp_Desai/'];
path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
numGtitle2 = f_ProxString_Curva(fid);
nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nGraf = nGraf{1};

% Numero steps provenientes del codigo
numSteps2 = f_ProxString_Curva(fid);
nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nSteps = nSteps{1}; 

% Lectura del t�tulo para eje de abscisas
% seccion = f_ProxString_Curva(fid);
% ejeX = f_tipoDat(seccion,fid);
% vs = f_ProxString_Curva(fid);

% Lectura del t�tulo para eje de aordenadas
% seccion = f_ProxString_Curva(fid);
% ejeY = f_tipoDat(seccion,fid);

% Lectura de datos XY
format = '%q %q';
xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
xy = xy{1};
xy = str2double(xy(:,:));
xy(:,2) = xy(:,2)*0.3048;
% Grafico 1
figure(1);
hold on
if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
    semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'k-d');
    set(gca,'XAxisLocation','top');
else
    plot(xy(2:nSteps,1),xy(2:nSteps,2),'k-d');
end

% grid on
% xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% titulo = strcat(ejeX,'-',ejeY);
% title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% hold on
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
    
    main_example_path =[pwd '/Examples/Desai/Desai_log/Comp_Desai/'];
    path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
    file_open = fullfile(path_file,file);
    fid = fopen(file_open,'r');    
       
% N�mero steps provenientes del codigo    
    numSteps = f_ProxString_Curva(fid);
    nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
    nSteps = nSteps{1}; 
    
% Lectura del t�tulo para eje de abscisas
%     seccion = f_ProxString_Curva(fid);
%     ejeX = f_tipoDat(seccion,fid);
%     vs = f_ProxString_Curva(fid);
    
% Lectura de datos XY 
    format = '%q %q';
    xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
    xy = xy{1};
    xy = str2double(xy(:,:));
    switch iGraf
        case {3,5,7,9,11,13,15}
            xy(:,2) = xy(:,2)*0.3048;
        case {2,4,6,8,10,12,14,16}
            xy(:,2) = xy(:,2)*(47.8803)/1000;
    end
    % Grafico 1
    figure(iGraf)
    hold on
    if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
        semilogx(xy(:,1),xy(:,2),'k-d');
        set(gca,'XAxisLocation','top');
    else
        plot(xy(:,1),xy(:,2),'k-d');
    end
    
%     grid on
%     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
%     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
%     titulo = strcat(ejeX,'-',ejeY);
%     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');

%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

end

% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % % DIRECTORIO EN LA LINEA INFERIOR
% % % main_example_path =[pwd '/Examples/Comp_DNS_Original/DNS_Comp_01/'];
% % % path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
% % % file_open = fullfile(path_file,file);
% % % fid = fopen(file_open,'r');
% % % 
% % % % Numero de graficos provenientes del codigo
% % % numGtitle3 = f_ProxString_Curva(fid);
% % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nGraf = nGraf{1};
% % % 
% % % % Numero steps provenientes del codigo
% % % numSteps3 = f_ProxString_Curva(fid);
% % % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nSteps = nSteps{1}; 
% % % 
% % % % Lectura del t�tulo para eje de abscisas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeX = f_tipoDat(seccion,fid);
% % % % vs = f_ProxString_Curva(fid);
% % % 
% % % % Lectura del t�tulo para eje de aordenadas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeY = f_tipoDat(seccion,fid);
% % % 
% % % % Lectura de datos XY
% % % format = '%q %q';
% % % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % xy = xy{1};
% % % xy = str2double(xy(:,:));
% % % 
% % % % Grafico 1
% % % figure(10);
% % % hold on
% % % if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %     semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'m-x');
% % %     set(gca,'XAxisLocation','top');
% % % else
% % %     plot(xy(2:nSteps,1),xy(2:nSteps,2),'m-x');
% % % end
% % % 
% % % % grid on
% % % % xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % titulo = strcat(ejeX,'-',ejeY);
% % % % title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % hold on
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % %*******************************************************************************
% % % %* LOOP PARA CASOS QUE nGraf > 1                                               *
% % % %*******************************************************************************
% % % 
% % % for iGraf = 2:2
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % %Archivo "Macro" de gr�ficos a ser abierto
% % %     actGraf = num2str(iGraf,'%03d');
% % %     macroI = strcat('Macro',actGraf,'.dat');
% % %     
% % %     main_example_path =[pwd '/Examples/Comp_DNS_Original/DNS_Comp_01/'];
% % %     path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
% % %     file_open = fullfile(path_file,file);
% % %     fid = fopen(file_open,'r');    
% % %        
% % % % N�mero steps provenientes del codigo    
% % %     numSteps = f_ProxString_Curva(fid);
% % %     nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % %     nSteps = nSteps{1}; 
% % %     
% % % % Lectura del t�tulo para eje de abscisas
% % % %     seccion = f_ProxString_Curva(fid);
% % % %     ejeX = f_tipoDat(seccion,fid);
% % % %     vs = f_ProxString_Curva(fid);
% % %     
% % % % Lectura de datos XY 
% % %     format = '%q %q';
% % %     xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % %     xy = xy{1};
% % %     xy = str2double(xy(:,:));
% % %     % Grafico 1
% % %     figure(11)
% % %     hold on
% % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %         semilogx(xy(:,1),xy(:,2),'m-x');
% % %         set(gca,'XAxisLocation','top');
% % %     else
% % %         plot(xy(:,1),xy(:,2),'m-x');
% % %     end
% % %     
% % % %     grid on
% % % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     titulo = strcat(ejeX,'-',ejeY);
% % % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % 
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % end
% % % 
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % % DIRECTORIO EN LA LINEA INFERIOR
% % % main_example_path =[pwd '/Examples/Comp_DNS_Original/MS_Comp2/'];
% % % path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
% % % file_open = fullfile(path_file,file);
% % % fid = fopen(file_open,'r');
% % % 
% % % % Numero de graficos provenientes del codigo
% % % numGtitle4 = f_ProxString_Curva(fid);
% % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nGraf = nGraf{1};
% % % 
% % % % Numero steps provenientes del codigo
% % % numSteps4 = f_ProxString_Curva(fid);
% % % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nSteps = nSteps{1}; 
% % % 
% % % % Lectura del t�tulo para eje de abscisas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeX = f_tipoDat(seccion,fid);
% % % % vs = f_ProxString_Curva(fid);
% % % 
% % % % Lectura del t�tulo para eje de aordenadas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeY = f_tipoDat(seccion,fid);
% % % 
% % % % Lectura de datos XY
% % % format = '%q %q';
% % % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % xy = xy{1};
% % % xy = str2double(xy(:,:));
% % % 
% % % % Grafico 1
% % % figure(10);
% % % hold on
% % % if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %     semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'k-o');
% % %     set(gca,'XAxisLocation','top');
% % % else
% % %     plot(xy(2:nSteps,1),xy(2:nSteps,2),'k-o');
% % % end
% % % 
% % % % grid on
% % % % xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % titulo = strcat(ejeX,'-',ejeY);
% % % % title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % hold on
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % %*******************************************************************************
% % % %* LOOP PARA CASOS QUE nGraf > 1                                               *
% % % %*******************************************************************************
% % % 
% % % for iGraf = 2:nGraf
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % %Archivo "Macro" de gr�ficos a ser abierto
% % %     actGraf = num2str(iGraf,'%03d');
% % %     macroI = strcat('Macro',actGraf,'.dat');
% % %     
% % %     main_example_path =[pwd '/Examples/Comp_DNS_Original/MS_Comp2/'];
% % %     path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
% % %     file_open = fullfile(path_file,file);
% % %     fid = fopen(file_open,'r');    
% % %        
% % % % N�mero steps provenientes del codigo    
% % %     numSteps = f_ProxString_Curva(fid);
% % %     nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % %     nSteps = nSteps{1}; 
% % %     
% % % % Lectura del t�tulo para eje de abscisas
% % % %     seccion = f_ProxString_Curva(fid);
% % % %     ejeX = f_tipoDat(seccion,fid);
% % % %     vs = f_ProxString_Curva(fid);
% % %     
% % % % Lectura de datos XY 
% % %     format = '%q %q';
% % %     xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % %     xy = xy{1};
% % %     xy = str2double(xy(:,:));
% % %     % Grafico 1
% % %     figure(11)
% % %     hold on
% % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %         semilogx(xy(:,1),xy(:,2),'k-o');
% % %         set(gca,'XAxisLocation','top');
% % %     else
% % %         plot(xy(:,1),xy(:,2),'k-o');
% % %     end
% % %     
% % % %     grid on
% % % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     titulo = strcat(ejeX,'-',ejeY);
% % % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % 
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % end
% % % %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % % %%
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % % DIRECTORIO EN LA LINEA INFERIOR
% % % main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE8_MatE4Imp/'];
% % % path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
% % % file_open = fullfile(path_file,file);
% % % fid = fopen(file_open,'r');
% % % 
% % % % Numero de graficos provenientes del codigo
% % % numGtitle5 = f_ProxString_Curva(fid);
% % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nGraf = nGraf{1};
% % % 
% % % % Numero steps provenientes del codigo
% % % numSteps5 = f_ProxString_Curva(fid);
% % % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nSteps = nSteps{1}; 
% % % 
% % % % Lectura del t�tulo para eje de abscisas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeX = f_tipoDat(seccion,fid);
% % % % vs = f_ProxString_Curva(fid);
% % % 
% % % % Lectura del t�tulo para eje de aordenadas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeY = f_tipoDat(seccion,fid);
% % % 
% % % % Lectura de datos XY
% % % format = '%q %q';
% % % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % xy = xy{1};
% % % xy = str2double(xy(:,:));
% % % 
% % % % Grafico 1
% % % figure(1);
% % % hold on
% % % if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %     semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'g-');
% % %     set(gca,'XAxisLocation','top');
% % % else
% % %     plot(xy(2:nSteps,1),xy(2:nSteps,2),'g-');
% % % end
% % % 
% % % % grid on
% % % % xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % titulo = strcat(ejeX,'-',ejeY);
% % % % title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % hold on
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % %*******************************************************************************
% % % %* LOOP PARA CASOS QUE nGraf > 1                                               *
% % % %*******************************************************************************
% % % 
% % % for iGraf = 2:nGraf
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % %Archivo "Macro" de gr�ficos a ser abierto
% % %     actGraf = num2str(iGraf,'%03d');
% % %     macroI = strcat('Macro',actGraf,'.dat');
% % %     
% % %     main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE8_MatE4Imp/'];
% % %     path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
% % %     file_open = fullfile(path_file,file);
% % %     fid = fopen(file_open,'r');    
% % %        
% % % % N�mero steps provenientes del codigo    
% % %     numSteps = f_ProxString_Curva(fid);
% % %     nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % %     nSteps = nSteps{1}; 
% % %     
% % % % Lectura del t�tulo para eje de abscisas
% % % %     seccion = f_ProxString_Curva(fid);
% % % %     ejeX = f_tipoDat(seccion,fid);
% % % %     vs = f_ProxString_Curva(fid);
% % %     
% % % % Lectura de datos XY 
% % %     format = '%q %q';
% % %     xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % %     xy = xy{1};
% % %     xy = str2double(xy(:,:));
% % %     % Grafico 1
% % %     figure(iGraf)
% % %     hold on
% % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %         semilogx(xy(:,1),xy(:,2),'g-');
% % %         set(gca,'XAxisLocation','top');
% % %     else
% % %         plot(xy(:,1),xy(:,2),'g-');
% % %     end
% % %     
% % % %     grid on
% % % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     titulo = strcat(ejeX,'-',ejeY);
% % % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % 
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % end
% % % %%
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % % DIRECTORIO EN LA LINEA INFERIOR
% % % main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE8_MatE16Imp/'];
% % % path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
% % % file_open = fullfile(path_file,file);
% % % fid = fopen(file_open,'r');
% % % 
% % % % Numero de graficos provenientes del codigo
% % % numGtitle6 = f_ProxString_Curva(fid);
% % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nGraf = nGraf{1};
% % % 
% % % % Numero steps provenientes del codigo
% % % numSteps6 = f_ProxString_Curva(fid);
% % % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nSteps = nSteps{1}; 
% % % 
% % % % Lectura del t�tulo para eje de abscisas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeX = f_tipoDat(seccion,fid);
% % % % vs = f_ProxString_Curva(fid);
% % % 
% % % % Lectura del t�tulo para eje de aordenadas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeY = f_tipoDat(seccion,fid);
% % % 
% % % % Lectura de datos XY
% % % format = '%q %q';
% % % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % xy = xy{1};
% % % xy = str2double(xy(:,:));
% % % 
% % % % Grafico 1
% % % figure(1);
% % % hold on
% % % if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %     semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'m-');
% % %     set(gca,'XAxisLocation','top');
% % % else
% % %     plot(xy(2:nSteps,1),xy(2:nSteps,2),'m-');
% % % end
% % % 
% % % % grid on
% % % % xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % titulo = strcat(ejeX,'-',ejeY);
% % % % title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % hold on
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % %*******************************************************************************
% % % %* LOOP PARA CASOS QUE nGraf > 1                                               *
% % % %*******************************************************************************
% % % 
% % % for iGraf = 2:nGraf
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % %Archivo "Macro" de gr�ficos a ser abierto
% % %     actGraf = num2str(iGraf,'%03d');
% % %     macroI = strcat('Macro',actGraf,'.dat');
% % %     
% % %     main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE8_MatE16Imp/'];
% % %     path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
% % %     file_open = fullfile(path_file,file);
% % %     fid = fopen(file_open,'r');    
% % %        
% % % % N�mero steps provenientes del codigo    
% % %     numSteps = f_ProxString_Curva(fid);
% % %     nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % %     nSteps = nSteps{1}; 
% % %     
% % % % Lectura del t�tulo para eje de abscisas
% % % %     seccion = f_ProxString_Curva(fid);
% % % %     ejeX = f_tipoDat(seccion,fid);
% % % %     vs = f_ProxString_Curva(fid);
% % %     
% % % % Lectura de datos XY 
% % %     format = '%q %q';
% % %     xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % %     xy = xy{1};
% % %     xy = str2double(xy(:,:));
% % %     % Grafico 1
% % %     figure(iGraf)
% % %     hold on
% % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %         semilogx(xy(:,1),xy(:,2),'m-');
% % %         set(gca,'XAxisLocation','top');
% % %     else
% % %         plot(xy(:,1),xy(:,2),'m-');
% % %     end
% % %     
% % % %     grid on
% % % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     titulo = strcat(ejeX,'-',ejeY);
% % % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % 
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % end
% % % %%
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % % DIRECTORIO EN LA LINEA INFERIOR
% % % main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE8_MatE36Imp/'];
% % % path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
% % % file_open = fullfile(path_file,file);
% % % fid = fopen(file_open,'r');
% % % 
% % % % Numero de graficos provenientes del codigo
% % % numGtitle7 = f_ProxString_Curva(fid);
% % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nGraf = nGraf{1};
% % % 
% % % % Numero steps provenientes del codigo
% % % numSteps7 = f_ProxString_Curva(fid);
% % % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nSteps = nSteps{1}; 
% % % 
% % % % Lectura del t�tulo para eje de abscisas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeX = f_tipoDat(seccion,fid);
% % % % vs = f_ProxString_Curva(fid);
% % % 
% % % % Lectura del t�tulo para eje de aordenadas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeY = f_tipoDat(seccion,fid);
% % % 
% % % % Lectura de datos XY
% % % format = '%q %q';
% % % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % xy = xy{1};
% % % xy = str2double(xy(:,:));
% % % 
% % % % Grafico 1
% % % figure(1);
% % % hold on
% % % if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %     semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'k-o');
% % %     set(gca,'XAxisLocation','top');
% % % else
% % %     plot(xy(2:nSteps,1),xy(2:nSteps,2),'k-o');
% % % end
% % % 
% % % % grid on
% % % % xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % titulo = strcat(ejeX,'-',ejeY);
% % % % title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % hold on
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % %*******************************************************************************
% % % %* LOOP PARA CASOS QUE nGraf > 1                                               *
% % % %*******************************************************************************
% % % 
% % % for iGraf = 2:nGraf
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % %Archivo "Macro" de gr�ficos a ser abierto
% % %     actGraf = num2str(iGraf,'%03d');
% % %     macroI = strcat('Macro',actGraf,'.dat');
% % %     
% % %     main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE8_MatE36Imp/'];
% % %     path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
% % %     file_open = fullfile(path_file,file);
% % %     fid = fopen(file_open,'r');    
% % %        
% % % % N�mero steps provenientes del codigo    
% % %     numSteps = f_ProxString_Curva(fid);
% % %     nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % %     nSteps = nSteps{1}; 
% % %     
% % % % Lectura del t�tulo para eje de abscisas
% % % %     seccion = f_ProxString_Curva(fid);
% % % %     ejeX = f_tipoDat(seccion,fid);
% % % %     vs = f_ProxString_Curva(fid);
% % %     
% % % % Lectura de datos XY 
% % %     format = '%q %q';
% % %     xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % %     xy = xy{1};
% % %     xy = str2double(xy(:,:));
% % %     % Grafico 1
% % %     figure(iGraf)
% % %     hold on
% % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %         semilogx(xy(:,1),xy(:,2),'k-o');
% % %         set(gca,'XAxisLocation','top');
% % %     else
% % %         plot(xy(:,1),xy(:,2),'k-o');
% % %     end
% % %     
% % % %     grid on
% % % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     titulo = strcat(ejeX,'-',ejeY);
% % % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % 
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % end
% % % % %%
% % % % %*******************************************************************************
% % % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % % %*******************************************************************************
% % % % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % % % DIRECTORIO EN LA LINEA INFERIOR
% % % % main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE16_2/'];
% % % % path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
% % % % file_open = fullfile(path_file,file);
% % % % fid = fopen(file_open,'r');
% % % % 
% % % % % Numero de graficos provenientes del codigo
% % % % numGtitle8 = f_ProxString_Curva(fid);
% % % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % % nGraf = nGraf{1};
% % % % 
% % % % % Numero steps provenientes del codigo
% % % % numSteps8 = f_ProxString_Curva(fid);
% % % % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % % nSteps = nSteps{1}; 
% % % % 
% % % % % Lectura del t�tulo para eje de abscisas
% % % % % seccion = f_ProxString_Curva(fid);
% % % % % ejeX = f_tipoDat(seccion,fid);
% % % % % vs = f_ProxString_Curva(fid);
% % % % 
% % % % % Lectura del t�tulo para eje de aordenadas
% % % % % seccion = f_ProxString_Curva(fid);
% % % % % ejeY = f_tipoDat(seccion,fid);
% % % % 
% % % % % Lectura de datos XY
% % % % format = '%q %q';
% % % % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % % xy = xy{1};
% % % % xy = str2double(xy(:,:));
% % % % 
% % % % % Grafico 1
% % % % figure(1);
% % % % hold on
% % % % if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % % %     semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'y-o');
% % % %     set(gca,'XAxisLocation','top');
% % % % else
% % % %     plot(xy(2:nSteps,1),xy(2:nSteps,2),'y-o');
% % % % end
% % % % 
% % % % % grid on
% % % % % xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % % ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % % titulo = strcat(ejeX,'-',ejeY);
% % % % % title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % % hold on
% % % % %*******************************************************************************
% % % % %* SE CIERRA ARCHIVO                                                           *
% % % % %*******************************************************************************
% % % % fclose(fid);
% % % % 
% % % % %*******************************************************************************
% % % % %* LOOP PARA CASOS QUE nGraf > 1                                               *
% % % % %*******************************************************************************
% % % % 
% % % % for iGraf = 2:nGraf
% % % % %*******************************************************************************
% % % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % % %*******************************************************************************
% % % % %Archivo "Macro" de gr�ficos a ser abierto
% % % %     actGraf = num2str(iGraf,'%03d');
% % % %     macroI = strcat('Macro',actGraf,'.dat');
% % % %     
% % % %     main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE16_2/'];
% % % %     path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
% % % %     file_open = fullfile(path_file,file);
% % % %     fid = fopen(file_open,'r');    
% % % %        
% % % % % N�mero steps provenientes del codigo    
% % % %     numSteps = f_ProxString_Curva(fid);
% % % %     nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % %     nSteps = nSteps{1}; 
% % % %     
% % % % % Lectura del t�tulo para eje de abscisas
% % % % %     seccion = f_ProxString_Curva(fid);
% % % % %     ejeX = f_tipoDat(seccion,fid);
% % % % %     vs = f_ProxString_Curva(fid);
% % % %     
% % % % % Lectura de datos XY 
% % % %     format = '%q %q';
% % % %     xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % %     xy = xy{1};
% % % %     xy = str2double(xy(:,:));
% % % %     % Grafico 1
% % % %     figure(iGraf)
% % % %     hold on
% % % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % % %         semilogx(xy(:,1),xy(:,2),'y-o');
% % % %         set(gca,'XAxisLocation','top');
% % % %     else
% % % %         plot(xy(:,1),xy(:,2),'y-o');
% % % %     end
% % % %     
% % % % %     grid on
% % % % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % %     titulo = strcat(ejeX,'-',ejeY);
% % % % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % 
% % % % %*******************************************************************************
% % % % %* SE CIERRA ARCHIVO                                                           *
% % % % %*******************************************************************************
% % % % fclose(fid);
% % % % 
% % % % end
% % % %%
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % % DIRECTORIO EN LA LINEA INFERIOR
% % % main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE16_3/'];
% % % path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
% % % file_open = fullfile(path_file,file);
% % % fid = fopen(file_open,'r');
% % % 
% % % % Numero de graficos provenientes del codigo
% % % numGtitle9 = f_ProxString_Curva(fid);
% % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nGraf = nGraf{1};
% % % 
% % % % Numero steps provenientes del codigo
% % % numSteps9 = f_ProxString_Curva(fid);
% % % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nSteps = nSteps{1}; 
% % % 
% % % % Lectura del t�tulo para eje de abscisas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeX = f_tipoDat(seccion,fid);
% % % % vs = f_ProxString_Curva(fid);
% % % 
% % % % Lectura del t�tulo para eje de aordenadas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeY = f_tipoDat(seccion,fid);
% % % 
% % % % Lectura de datos XY
% % % format = '%q %q';
% % % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % xy = xy{1};
% % % xy = str2double(xy(:,:));
% % % 
% % % % Grafico 1
% % % figure(1);
% % % hold on
% % % if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %     semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'r-o');
% % %     set(gca,'XAxisLocation','top');
% % % else
% % %     plot(xy(2:nSteps,1),xy(2:nSteps,2),'r-o');
% % % end
% % % 
% % % % grid on
% % % % xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % titulo = strcat(ejeX,'-',ejeY);
% % % % title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % hold on
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % %*******************************************************************************
% % % %* LOOP PARA CASOS QUE nGraf > 1                                               *
% % % %*******************************************************************************
% % % 
% % % for iGraf = 2:nGraf
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % %Archivo "Macro" de gr�ficos a ser abierto
% % %     actGraf = num2str(iGraf,'%03d');
% % %     macroI = strcat('Macro',actGraf,'.dat');
% % %     
% % %     main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE16_3/'];
% % %     path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
% % %     file_open = fullfile(path_file,file);
% % %     fid = fopen(file_open,'r');    
% % %        
% % % % N�mero steps provenientes del codigo    
% % %     numSteps = f_ProxString_Curva(fid);
% % %     nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % %     nSteps = nSteps{1}; 
% % %     
% % % % Lectura del t�tulo para eje de abscisas
% % % %     seccion = f_ProxString_Curva(fid);
% % % %     ejeX = f_tipoDat(seccion,fid);
% % % %     vs = f_ProxString_Curva(fid);
% % %     
% % % % Lectura de datos XY 
% % %     format = '%q %q';
% % %     xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % %     xy = xy{1};
% % %     xy = str2double(xy(:,:));
% % %     % Grafico 1
% % %     figure(iGraf)
% % %     hold on
% % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %         semilogx(xy(:,1),xy(:,2),'r-o');
% % %         set(gca,'XAxisLocation','top');
% % %     else
% % %         plot(xy(:,1),xy(:,2),'r-o');
% % %     end
% % %     
% % % %     grid on
% % % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     titulo = strcat(ejeX,'-',ejeY);
% % % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % 
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % end
% % % %%
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % % PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% % % % DIRECTORIO EN LA LINEA INFERIOR
% % % main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE16_4/'];
% % % path_file= main_example_path; file = 'Macro001.dat' ; %'\Examples\MS_Comp\'
% % % file_open = fullfile(path_file,file);
% % % fid = fopen(file_open,'r');
% % % 
% % % % Numero de graficos provenientes del codigo
% % % numGtitle10 = f_ProxString_Curva(fid);
% % % nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nGraf = nGraf{1};
% % % 
% % % % Numero steps provenientes del codigo
% % % numSteps10 = f_ProxString_Curva(fid);
% % % nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % % nSteps = nSteps{1}; 
% % % 
% % % % Lectura del t�tulo para eje de abscisas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeX = f_tipoDat(seccion,fid);
% % % % vs = f_ProxString_Curva(fid);
% % % 
% % % % Lectura del t�tulo para eje de aordenadas
% % % % seccion = f_ProxString_Curva(fid);
% % % % ejeY = f_tipoDat(seccion,fid);
% % % 
% % % % Lectura de datos XY
% % % format = '%q %q';
% % % xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % % xy = xy{1};
% % % xy = str2double(xy(:,:));
% % % 
% % % % Grafico 1
% % % figure(1);
% % % hold on
% % % if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %     semilogx(xy(2:nSteps,1),xy(2:nSteps,2),'r-o');
% % %     set(gca,'XAxisLocation','top');
% % % else
% % %     plot(xy(2:nSteps,1),xy(2:nSteps,2),'r-o');
% % % end
% % % 
% % % % grid on
% % % % xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % % ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % % titulo = strcat(ejeX,'-',ejeY);
% % % % title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % % hold on
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % %*******************************************************************************
% % % %* LOOP PARA CASOS QUE nGraf > 1                                               *
% % % %*******************************************************************************
% % % 
% % % for iGraf = 2:nGraf
% % % %*******************************************************************************
% % % %* APERTURA DEL ARCHIVO DE DATOS                                               *
% % % %*******************************************************************************
% % % %Archivo "Macro" de gr�ficos a ser abierto
% % %     actGraf = num2str(iGraf,'%03d');
% % %     macroI = strcat('Macro',actGraf,'.dat');
% % %     
% % %     main_example_path =[pwd '/Examples/Test04_05/Tlog2Mat/MS_RVE16_4/'];
% % %     path_file= main_example_path; file = macroI ; %'\Examples\MS_Comp\'
% % %     file_open = fullfile(path_file,file);
% % %     fid = fopen(file_open,'r');    
% % %        
% % % % N�mero steps provenientes del codigo    
% % %     numSteps = f_ProxString_Curva(fid);
% % %     nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% % %     nSteps = nSteps{1}; 
% % %     
% % % % Lectura del t�tulo para eje de abscisas
% % % %     seccion = f_ProxString_Curva(fid);
% % % %     ejeX = f_tipoDat(seccion,fid);
% % % %     vs = f_ProxString_Curva(fid);
% % %     
% % % % Lectura de datos XY 
% % %     format = '%q %q';
% % %     xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
% % %     xy = xy{1};
% % %     xy = str2double(xy(:,:));
% % %     % Grafico 1
% % %     figure(iGraf)
% % %     hold on
% % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %         semilogx(xy(:,1),xy(:,2),'r-o');
% % %         set(gca,'XAxisLocation','top');
% % %     else
% % %         plot(xy(:,1),xy(:,2),'r-o');
% % %     end
% % %     
% % % %     grid on
% % % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % % %     titulo = strcat(ejeX,'-',ejeY);
% % % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');
% % % 
% % % %*******************************************************************************
% % % %* SE CIERRA ARCHIVO                                                           *
% % % %*******************************************************************************
% % % fclose(fid);
% % % 
% % % end