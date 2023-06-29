% Este script permite graficar las curvas dadas en los archivos de salida
% de codigo Multiscale_Application

clear;
clc;
close all
%*******************************************************************************
%* RVE 1x1 - FOE                                                               *
%*******************************************************************************
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Periodico/Terzaghi/Monoescala_Ref2/'];
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
z_Hn2 = [0 (10-9.9775)/10,... 
    (10-9.9550)/10,...  
    (10-9.9121)/10,...  
    (10-9.8692)/10,...  
    (10-9.7911)/10,...  
    (10-9.7130)/10,...  
    (10-9.5782)/10,...  
    (10-9.4434)/10,...  
    (10-9.2227)/10,...  
    (10-9.0021)/10,...  
    (10-8.3202)/10,...  
    (10-7.3282)/10,...  
    (10-5.9756)/10,...  
    (10-4.2540)/10,...  
    (10-2.2195)/10,... 
    1];
st2 = size(z_Hn2,2);
p_p0_Tv1 = zeros(st2,1);
p_p0_Tv2 = zeros(st2,1);
p_p0_Tv3 = zeros(st2,1);
inte2=1/st2;

for iGraf = 2:nGraf
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
%Archivo "Macro" de graficos a ser abierto
    actGraf = num2str(iGraf,'%03d');
    macroI = strcat('Macro',actGraf,'.dat');

    main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Periodico/Terzaghi/Monoescala_Ref2'];
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
    
    p_p0_Tv1(iGraf-1,1) = xy(11,2);
    p_p0_Tv2(iGraf-1,1) = xy(15,2);
    p_p0_Tv3(iGraf-1,1) = xy(20,2);
 
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

end
%*******************************************************************************
%* RVE 1X1 SOE                                             *
%*******************************************************************************
%################################################################################
%################################################################################
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Periodico/Terzaghi/Monoescala_Ref2'];
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
z_Hn = [0 (10-9.797)/10 (10-9.492)/10 (10-9.036)/10 (10-8.351)/10 (10-7.323)/10 (10-5.781)/10 (10-3.469)/10 1];
% k2n = size(z_Hn,2);
st = size(z_Hn,2);
p_p0_Tv1_2 = zeros(st2,1);
p_p0_Tv2_2 = zeros(st2,1);
p_p0_Tv3_2= zeros(st2,1);
inte=1/st; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&&&&&&
for iGraf = 2:nGraf
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
%Archivo "Macro" de graficos a ser abierto
    actGraf = num2str(iGraf,'%03d');
    macroI = strcat('Macro',actGraf,'.dat');

    main_example_path =['/home/reiadrian/workspace/Multiscale_PorousMedia/Ejemplos/Homogeneos/Periodico/Terzaghi/Monoescala_Ref2'];
    path_file= main_example_path; file = macroI ; 
    file_open = fullfile(path_file,file);
    fid = fopen(file_open,'r');    
       
% Numero steps provenientes del codigo    
    numSteps = f_ProxString_Curva(fid);
    nSteps = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
%     nSteps = nSteps{1}; 
    nSteps =1179; 
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
    
    p_p0_Tv1_2(iGraf-1,1) = xy(11,2);
    p_p0_Tv2_2(iGraf-1,1) = xy(15,2);
    p_p0_Tv3_2(iGraf-1,1) = xy(20,2);

%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

end
%################################################################################
%################################################################################
%Calculo de cv 
gamma_w=0.981;%kN/m3
E=100;%kPa
nu=0.3;
e0=0.6;
k=0.0000864;%m/dia
factor=1-((2*(nu^2))/(1-nu));
Em = E/factor;
K=(3*E)/(3*(1-2*nu));
mv=(1/Em);%/(1+e0);
cv=k/(mv*gamma_w);
H = 10;
cv_H2 = cv/(H^2);
Tv = [cv_H2*1 cv_H2*5 cv_H2*10] ;
% Tv = [cv_H2*10] ;
z_H = 0:0.01:1;
m = 0:100000;
k1 = size(Tv,2);
k2 = size(z_H,2);
k3 = size(m,2);
p_p0 = zeros(k2,1,k1);
for ik1=1:k1
    for ik2=1:k2
        for ik3=1:k3
            M = 0.5*pi()*(2*m(ik3)+1);
            Mz_H = M*z_H(ik2);
            sMz_H = sin(Mz_H);
            eTv = exp(1)^(-(M^2)*Tv(ik1));
            p_p0(ik2,1,ik1) = p_p0(ik2,1,ik1) + ((2/M)*sMz_H*eTv);
        end
    end

end
pos = [1 2 3 4 5 7 14 15 17];
figure(1)

h = gobjects(3,2); 
h(1,1) =plot(p_p0_Tv1,z_Hn2,'r-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
axis([0 1.04 -0.099 1])
hold on
h(2,1)=plot(p_p0_Tv1_2,z_Hn2,'b-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
hold on
h(3,1) = plot(p_p0(:,1,1),z_H,'k-','LineWidth',4,'MarkerSize',10,'DisplayName', 'Closed-form');
% hold on
% h(4,1) =plot(p_p0_Tv1(pos),z_Hn2(pos),'r*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
% hold on
% h(5,1)=plot(p_p0_Tv1_2(pos),z_Hn2(pos),'ks','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
% hold on

grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',24,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
lg  = legend(h(1:3,1)); 
lg.FontSize = 24;

figure(2)
h(1,2) = plot(p_p0_Tv1,z_Hn2,'r-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
axis([0.8 1.04 0 0.4])
hold on
h(2,2) =plot(p_p0_Tv1_2,z_Hn2,'b-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
hold on
h(3,2) = plot(p_p0(:,1,1),z_H,'k-','LineWidth',4,'MarkerSize',10,'DisplayName', 'Closed-form');
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',24,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;

%################################################################################################
figure(3)
h2 = gobjects(3,2); 
% subplot(2,1,1);
h2(1,1) = plot(p_p0_Tv2,z_Hn2,'r-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
axis([0 1.04 -0.099 1])
hold on
h2(2,1) =plot(p_p0_Tv2_2,z_Hn2,'b-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
hold on
h2(3,1) = plot(p_p0(:,1,2),z_H,'k-','LineWidth',4,'MarkerSize',10,'DisplayName', 'Closed-form');
% hold on
% h2(4,1) =plot(p_p0_Tv2(pos),z_Hn2(pos),'r*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
% hold on
% h2(5,1)=plot(p_p0_Tv2_2(pos),z_Hn2(pos),'ks','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
% hold on
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',24,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
lg2  = legend(h2(1:3,1)); 
lg2.FontSize = 24;

figure(4)
h2(1,2) = plot(p_p0_Tv2,z_Hn2,'r-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
hold on
h2(2,2) = plot(p_p0_Tv2_2,z_Hn2,'b-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
hold on
h2(3,2) = plot(p_p0(:,1,2),z_H,'k-','LineWidth',4,'MarkerSize',10,'DisplayName', 'Closed-form');
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',24,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
%################################################################################################
figure(5)
h3 = gobjects(3,2);
h3(1,1) = plot(p_p0_Tv3,z_Hn2,'r-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
axis([0 1.04 -0.099 1])
hold on
h3(2,1) =plot(p_p0_Tv3_2,z_Hn2,'b-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
hold on
h3(3,1) =plot(p_p0(:,1,3),z_H,'k-','LineWidth',4,'MarkerSize',10,'DisplayName', 'Closed-form');
% hold on
% h2(4,1) =plot(p_p0_Tv3(pos),z_Hn2(pos),'r*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
% hold on
% h2(5,1)=plot(p_p0_Tv3_2(pos),z_Hn2(pos),'ks','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
% hold on
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',24,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;
lg3  = legend(h3(1:3,1)); 
lg3.FontSize = 24;

figure(6)
h3(1,2) = plot(p_p0_Tv3,z_Hn2,'r-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'RCE');
axis([0.8 1.04 0 0.4])
hold on
h3(2,2) = plot(p_p0_Tv3_2,z_Hn2,'b-*','LineWidth',4,'MarkerSize',10,'DisplayName', 'OCE');
hold on
h3(3,2) = plot(p_p0(:,1,3),z_H,'k-','LineWidth',4,'MarkerSize',10,'DisplayName', 'Closed-form');
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',24,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',24,'FontWeight','bold','Color','k');
ax.FontSize = 24;
ay.FontSize = 24;

