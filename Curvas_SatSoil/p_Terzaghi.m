% % % % %Calculo de poro-presiones con la formula de Terzaghi para diferentes
% % % % %tiempos (y una misma profundidad)
% % % % clear all
% % % % close all
% % % % clc
% % % % %Calculo de cv 
% % % % gamma_w=9.81;%kN/m3
% % % % E=100;%kPa
% % % % nu=0.3;
% % % % e0=0.6;
% % % % k=0.0000864;%m/dia
% % % % factor=1-((2*(nu^2))/(1-nu));
% % % % Em = E/factor;
% % % % K=(3*E)/(3*(1-2*nu));
% % % % mv=(1/Em);%/(1+e0);
% % % % cvd=k/(mv*gamma_w);
% % % % 
% % % % % cv_H2 = 0.10496e-4;
% % % % 
% % % % % M = 0.5*pi()*(2*m+1);
% % % % z = 10-0.964;
% % % % H = 10;
% % % % p0 = gamma_w*z;
% % % % % z = H-z1;
% % % % cv_H2 = cvd/H^2;
% % % % z_H = z/H;
% % % % t = 0:0.1:60000;
% % % % t =t';
% % % % m = 1:1000;
% % % % m=m';
% % % % k2 = size(t,1);
% % % % k3 = size(m,1);
% % % % p_p0 = zeros(k2,1);
% % % % for ik2=1:k2
% % % %     Tv = cv_H2*t(ik2);
% % % %     for ik3=1:k3
% % % %         M = 0.5*pi()*(2*m(ik3)+1);
% % % %         Mz_H = M*z_H;
% % % %         sMz_H = sin(Mz_H);       
% % % %         eTv = exp(1)^(-(M^2)*Tv);
% % % %         p_p0(ik2,1) = p_p0(ik2,1) + ((2*p0/M)*sMz_H*eTv);
% % % %     end
% % % % end
% % % % figure(2)
% % % % semilogx(t,p_p0)


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
main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Rel_Const_Actual/RVE1x1_UCR/'];
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
p_p0_Tv1 = zeros(st,1);
p_p0_Tv2 = zeros(st,1);
p_p0_Tv3 = zeros(st,1);
inte=1/st; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&&&&&&
for iGraf = 2:2:19
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
%Archivo "Macro" de graficos a ser abierto
    actGraf = num2str(iGraf,'%03d');
    macroI = strcat('Macro',actGraf,'.dat');

    main_example_path =['/home/reiadrian/workspace/LAGRANGE/Ejemplos/Rel_Const_Actual/RVE1x1_UCR/'];
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
    
    p_p0_Tv1(iGraf/2,1) = xy(21,2);
    p_p0_Tv2(iGraf/2,1) = xy(51,2);
    p_p0_Tv3(iGraf/2,1) = xy(101,2);
 
%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

end
%################################################################################
%################################################################################
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
% PARA MODIFIFICAR EL ARCHIVO DE INGRESO SOLO HACE FALTA CAMBIAR EL
% DIRECTORIO EN LA LINEA INFERIOR
main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Rel_Const_Anterior/RVE1x1_OCR/'];
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
p_p0_Tv1_2 = zeros(st,1);
p_p0_Tv2_2 = zeros(st,1);
p_p0_Tv3_2= zeros(st,1);
inte=1/st; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!&&&&&&&&&&&&&&&&&&&&&&&
for iGraf = 2:2:19
%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
%Archivo "Macro" de graficos a ser abierto
    actGraf = num2str(iGraf,'%03d');
    macroI = strcat('Macro',actGraf,'.dat');

    main_example_path =['/home/reiadrian/workspace/LAGRANGE2/Ejemplos/Rel_Const_Anterior/RVE1x1_OCR/'];
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
    
    p_p0_Tv1_2(iGraf/2,1) = xy(21,2);
    p_p0_Tv2_2(iGraf/2,1) = xy(51,2);
    p_p0_Tv3_2(iGraf/2,1) = xy(101,2);

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
Tv = [cv_H2*2 cv_H2*5 cv_H2*10] ;
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

figure(1)
h = gobjects(3,2); 

subplot(2,1,1);
h(1,1) =plot(p_p0_Tv1,z_Hn,'r-*','LineWidth',2,'MarkerSize',10,'DisplayName', 'UCR');
axis([0 1.04 0 1])
hold on
h(2,1)=plot(p_p0_Tv1_2,z_Hn,'k-s','LineWidth',2,'MarkerSize',10,'DisplayName', 'OCR');
axis([0 1.04 0 1])
hold on
h(3,1) = plot(p_p0(:,1,1),z_H,'b','LineWidth',2,'MarkerSize',10,'DisplayName', 'Closed-form');
axis([0 1.04 0 1])
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',12,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',12,'FontWeight','bold','Color','k');
ax.FontSize = 12;
ay.FontSize = 12;
% title('T_{v} =  2.37e-04 - t = 2s','FontSize',12)

subplot(2,1,2); 
h(1,2) = plot(p_p0_Tv1,z_Hn,'r-*','LineWidth',2,'MarkerSize',10,'DisplayName', 'UCR');
axis([0.8 1.04 0 0.4])
hold on
h(2,2) =plot(p_p0_Tv1_2,z_Hn,'k-s','LineWidth',2,'MarkerSize',10,'DisplayName', 'OCR');
axis([0.8 1.04 0 0.4])
hold on
h(3,2) = plot(p_p0(:,1,1),z_H,'b','LineWidth',2,'MarkerSize',10,'DisplayName', 'Closed-form');
axis([0.8 1.04 0 0.4])
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',12,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',12,'FontWeight','bold','Color','k');
ax.FontSize = 12;
ay.FontSize = 12;
% title('T_{v} =  2.37e-04 - t = 2s','FontSize',12)
lg  = legend(h(1:3,1)); 
lg.Position(1:2) = [.47 0.005];
lg.FontSize = 14;
%################################################################################################
figure(2)
h2 = gobjects(3,2); 
subplot(2,1,1);
h2(1,1) = plot(p_p0_Tv2,z_Hn,'r-*','LineWidth',2,'MarkerSize',10,'DisplayName', 'UCR');
axis([0 1.02 0 1])
hold on
h2(2,1) =plot(p_p0_Tv2_2,z_Hn,'k-s','LineWidth',2,'MarkerSize',10,'DisplayName', 'OCR');
axis([0 1.02 0 1])
hold on
h2(3,1) = plot(p_p0(:,1,2),z_H,'b','LineWidth',2,'MarkerSize',10,'DisplayName', 'Closed-form');
axis([0 1.02 0 1])
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',12,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',12,'FontWeight','bold','Color','k');
ax.FontSize = 12;
ay.FontSize = 12;
% title('T_{v} =  5.92e-04 - t = 5s','FontSize',12)

subplot(2,1,2); 
h2(1,2) = plot(p_p0_Tv2,z_Hn,'r-*','LineWidth',2,'MarkerSize',10,'DisplayName', 'UCR');
axis([0.8 1.02 0 0.4])
hold on
h2(2,2) = plot(p_p0_Tv2_2,z_Hn,'k-s','LineWidth',2,'MarkerSize',10,'DisplayName', 'OCR');
axis([0.8 1.02 0 0.4])
hold on
h2(3,2) = plot(p_p0(:,1,2),z_H,'b','LineWidth',2,'MarkerSize',10,'DisplayName', 'Closed-form');
axis([0.8 1.02 0 0.4])
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',12,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',12,'FontWeight','bold','Color','k');
ax.FontSize = 12;
ay.FontSize = 12;
%title('T_{v} =  5.92e-04 - t = 5s','FontSize',12)
lg2  = legend(h2(1:3,1)); 
lg2.Position(1:2) = [.47 0.005];
lg2.FontSize = 14;
%################################################################################################
figure(3)
h3 = gobjects(3,2);
subplot(2,1,1);
h3(1,1) = plot(p_p0_Tv3,z_Hn,'r-*','LineWidth',2,'MarkerSize',10,'DisplayName', 'UCR');
axis([0 1.01 0 1])
hold on
h3(2,1) =plot(p_p0_Tv3_2,z_Hn,'k-s','LineWidth',2,'MarkerSize',10,'DisplayName', 'OCR');
axis([0 1.01 0 1])
hold on
h3(3,1) =plot(p_p0(:,1,3),z_H,'b','LineWidth',2,'MarkerSize',10,'DisplayName', 'Closed-form');
axis([0 1.01 0 1])
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',12,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',12,'FontWeight','bold','Color','k');
ax.FontSize = 12;
ay.FontSize = 12;
% title('T_{v} =  0.12e-02 - t = 10s','FontSize',12)

subplot(2,1,2); 
h3(1,2) = plot(p_p0_Tv3,z_Hn,'r-*','LineWidth',2,'MarkerSize',10,'DisplayName', 'UCR');
axis([0.8 1.01 0 0.4])
hold on
h3(2,2) = plot(p_p0_Tv3_2,z_Hn,'k-s','LineWidth',2,'MarkerSize',10,'DisplayName', 'OCR');
axis([0.8 1.01 0 0.4])
hold on
h3(3,2) = plot(p_p0(:,1,3),z_H,'b','LineWidth',2,'MarkerSize',10,'DisplayName', 'Closed-form');
axis([0.8 1.01 0 0.4])
grid on
grid minor
ax = gca;
ax.YDir = 'reverse';
ylabel('z/H','FontSize',12,'FontWeight','bold','Color','k');
xlabel('p/p_{0}','FontSize',12,'FontWeight','bold','Color','k');
ax.FontSize = 12;
ay.FontSize = 12;
% title('T_{v} =  0.12e-02 - t = 10s','FontSize',12)
lg3  = legend(h3(1:3,1)); 
lg3.Position(1:2) = [.47 0.005];
lg3.FontSize = 14;

% fig.Position(4) = fig.Position(5) - 250;
% add legend
% lg  = legend(h(1:3,1)); 
% lg.Position(1:2) = [.47 0.005];
% lg.FontSize = 14;
% plot(p_p0_Tv2,z_Hn,'b-*')
% hold on
% 
% plot(p_p0_Tv3,z_Hn,'k-*')
% hold on
% 
% plot(p_p0_Tv1_2,z_Hn,'r-s')
% hold on
% 
% plot(p_p0_Tv2_2,z_Hn,'b-s')
% hold on
% 
% plot(p_p0_Tv3_2,z_Hn,'k-s')
% hold on


