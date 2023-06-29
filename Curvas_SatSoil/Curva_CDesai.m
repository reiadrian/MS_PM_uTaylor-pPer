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
main_example_path =[pwd '/Examples/Desai/Comp_Desai_SI/'];
path_file= main_example_path; file = 'Macro001.dat' ;
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Numero de graficos provenientes del codigo
numGtitle = f_ProxString_Curva(fid);
nGraf = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
nGraf = nGraf{1}; %El primer archivo de curvas tiene el numero total de curvas impresas

%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* LOOP PARA CASOS QUE nGraf > 1                                               *
%*******************************************************************************
z_Hn = [0 1/7 2/7 3/7 4/7 5/7 6/7 1];
k2n = size(z_Hn,2);
p_p0_Tv1 = zeros(k2n,1);
p_p0_Tv2 = zeros(k2n,1);
p_p0_Tv3 = zeros(k2n,1);
% p_p0_Tv4 = zeros(k2n,1);
% p_p0_Tv5 = zeros(k2n,1);
% p_p0_Tv6 = zeros(k2n,1);
% p_p0_Tv7 = zeros(k2n,1);
% p_p0_Tv8 = zeros(k2n,1);
% p_p0_Tv9 = zeros(k2n,1);
% p_p0_Tv10 = zeros(k2n,1);
for iGraf = 2:2:nGraf
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
%Archivo "Macro" de graficos a ser abierto
    actGraf = num2str(iGraf,'%03d');
    macroI = strcat('Macro',actGraf,'.dat');

    main_example_path =[pwd '/Examples/Desai/Comp_Desai_SI/'];
    path_file= main_example_path; file = macroI ; 
    file_open = fullfile(path_file,file);
    fid = fopen(file_open,'r');    
       
% Numero steps provenientes del codigo    
    numSteps = f_ProxString_Curva(fid);
    nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
    nSteps = nSteps{1}; 
    
% Lectura del titulo para eje de abscisas
    seccion = f_ProxString_Curva(fid);
    ejeX = f_tipoDat(seccion,fid);
    vs = f_ProxString_Curva(fid);
    
% Lectura del titulo para eje de aordenadas
    seccion = f_ProxString_Curva(fid);
    ejeY = f_tipoDat(seccion,fid);
    
% Lectura de datos XY 
    format = '%q %q';
    xy = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
    xy = xy{1};
    xy = str2double(xy(:,:));
%     p_p0_Tv1(iGraf/2,1) = xy(132,2);
%     p_p0_Tv2(iGraf/2,1) = xy(227,2);
%     p_p0_Tv3(iGraf/2,1) = xy(323,2);
%     p_p0_Tv4(iGraf/2,1) = xy(345,2);
%     p_p0_Tv5(iGraf/2,1) = xy(355,2);
%     p_p0_Tv6(iGraf/2,1) = xy(364,2);
%     p_p0_Tv7(iGraf/2,1) = xy(374,2);
%     p_p0_Tv8(iGraf/2,1) = xy(383,2);
%     p_p0_Tv9(iGraf/2,1) = xy(393,2);
%     p_p0_Tv10(iGraf/2,1) = xy(402,2);
%     if xy(2,2)==0
%         p_p0_Tv1(iGraf/2,1) = 0;
%         p_p0_Tv2(iGraf/2,1) = 0;
%         p_p0_Tv3(iGraf/2,1) = 0;
%     else
        p_p0_Tv1(iGraf/2,1) = xy(390,2)/(47.8803/1000);
%         p_p0_Tv2(iGraf/2,1) = xy(339,2)/(47.8803/1000);
%         p_p0_Tv3(iGraf/2,1) = xy(379,2)/(47.8803/1000);
%         p_p0_Tv1(iGraf/2,1) = xy(120,2);
%         p_p0_Tv2(iGraf/2,1) = xy(339,2);
%         p_p0_Tv3(iGraf/2,1) = xy(379,2);
%         p_p0_Tv1(iGraf/2,1) = xy(45,2)/xy(2,2);
%         p_p0_Tv2(iGraf/2,1) = xy(69,2)/xy(2,2);
%         p_p0_Tv3(iGraf/2,1) = xy(109,2)/xy(2,2);
%     end
    a= 1;
% % %     % Grafico 1
% % %     figure(iGraf)
% % %     if strcmp(ejeX,'logT') || strcmp(ejeY,'logT')
% % %         semilogx(xy(:,1),xy(:,2),'r-');
% % %         set(gca,'XAxisLocation','top');
% % %     else
% % %         plot(xy(:,1),xy(:,2),'r-');
% % %     end
% % %     
% % %     grid on
% % %     xlabel(ejeX,'FontSize',12,'FontWeight','bold','Color','k');
% % %     ylabel(ejeY,'FontSize',12,'FontWeight','bold','Color','k');
% % %     titulo = strcat(ejeX,'-',ejeY);
% % %     title(titulo,'FontSize',14,'FontWeight','bold','Color','k');

%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

end
Tv = [0.088 0.34 0.76];
% Tv = [0.088 0.1 0.2 0.3 0.34 0.4 0.5 0.6 0.7 0.76 0.8 0.9];
% Tv = 0.1:0.1:1;
z_H = 0:0.05:1;
m = 0:100000;
k1 = size(Tv,2);
k2 = size(z_H,2);
k3 = size(m,2);

for ik1=1:k1
    p_p0 = zeros(k2,1);
    for ik2=1:k2
        for ik3=1:k3
            M = 0.5*pi()*(2*m(ik3)+1);
            Mz_H = M*z_H(ik2);
            sMz_H = sin(Mz_H);
            eTv = exp(1)^(-(M^2)*Tv(ik1));
            p_p0(ik2,1) = p_p0(ik2,1) + ((2/M)*sMz_H*eTv);
        end
    end
    figure(1)
    plot(p_p0,z_H)
    hold on
end
plot(p_p0_Tv1,z_Hn,'k-*')
hold on
plot(p_p0_Tv2,z_Hn,'k-*')
hold on
plot(p_p0_Tv3,z_Hn,'k-*')
hold on
% plot(p_p0_Tv4,z_Hn,'k-*')
% hold on
% plot(p_p0_Tv5,z_Hn,'k-*')
% hold on
% plot(p_p0_Tv6,z_Hn,'k-*')
% hold on
% plot(p_p0_Tv7,z_Hn,'k-*')
% hold on
% plot(p_p0_Tv8,z_Hn,'k-*')
% hold on
% plot(p_p0_Tv9,z_Hn,'k-*')
% hold on
% plot(p_p0_Tv10,z_Hn,'k-*')
% hold on
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
