%% RVE_MetallicAlloy3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = [0.13     0   0  0  0; % 1
     0.078 0.078  0  0  0; % 2 
     0     0.13   0  0  0; % 3
    -0.13  0.13   0  0  0; % 4
    -0.13  0      0  0  0; % 5
    -0.095 -0.095 0  0  0; % 6
     0    -0.20   0  0  0; % 7
     0.4 -0.4   0  0  0; % 8
     0     0      0 -0.25 -0.25; % 9
     0.13  0      0 -0.13 -0.13; % 10
     0.07  0.07   0 -0.13 -0.13; % 11
     0     0.13   0 -0.13 -0.13; % 12
    -0.13  0.13   0 -0.13 -0.13;
    -0.13  0      0 -0.13 -0.13; 
    -0.09 -0.09   0 -0.09 -0.09; % 15
     0    -0.13   0 -0.13 -0.13; % 16
     0.13 -0.13   0 -0.13 -0.13;
     0     0      0  0.13  0.13; % 18
     0.13  0      0  0.13  0.13;
     0.06  0.06   0  0.13  0.13; % 20
     0     0.13   0  0.13  0.13;
    -0.13  0.13   0  0.13  0.13;
    -0.13  0      0  0.13  0.13;
    -0.10 -0.10   0  0.13  0.13;
     0    -0.13   0  0.13  0.13; % 25
     0.13 -0.13   0  0.13  0.13;
     0.0   0.0    0  0.20  0.0 ; 
     0     0      0  0    -0.20 ] ;
% plot3(t(:,1),t(:,2),t(:,4),'o') ; grid
for itra = 1: size(t,1)
    Macro0{itra}  = FACT*t(itra,:)' ;
end
%% RVE_GrapInclu03 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macro0{1}  = FACT*[ 0.097; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{2}  = FACT*[ 0.063; 0.063; 0.00; 0.00; 0.00] ; 
% Macro0{3}  = FACT*[ 0.00; 0.10; 0.00; 0.00; 0.00] ; 
% Macro0{4}  = FACT*[-0.25; 0.10; 0.00; 0.00; 0.00] ;
% Macro0{5}  = FACT*[-0.20; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{6}  = FACT*[-0.06;-0.10; 0.00; 0.00; 0.00] ;
% Macro0{7}  = FACT*[ 0.00;-0.15; 0.00; 0.00; 0.00] ;
% Macro0{8}  = FACT*[ 0.19;-0.10; 0.00; 0.00; 0.00] ;
% Macro0{9}  = FACT*[ 0.00; 0.00; 0.00; 0.30; 0.00] ;
% Macro0{10} = FACT*[ 0.14; 0.00; 0.00; 0.33; 0.00] ;
% Macro0{11} = 1.60*[ 0.04; 0.04; 0.00; 0.10; 0.00] ;
% Macro0{12} = 1.50*[ 0.00; 0.07; 0.00; 0.10; 0.00] ;
% Macro0{13} = FACT*[-0.17; 0.34; 0.00; 0.34; 0.00] ;
% Macro0{14} = FACT*[-0.22; 0.00; 0.00; 0.44; 0.00] ;
% Macro0{15} = FACT*[-0.10;-0.10; 0.00; 0.15; 0.00] ;
% Macro0{16} = 2.50*[ 0.00;-0.10; 0.00; 0.20; 0.00] ;
% Macro0{17} = FACT*[ 0.19;-0.10; 0.00; 0.35; 0.00] ;
%% RVE_GrapInclu02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macro0{1}  = FACT*[ 0.10; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{2}  = FACT*[ 0.066; 0.066; 0.00; 0.00; 0.00] ; 
% Macro0{3}  = FACT*[ 0.00; 0.10; 0.00; 0.00; 0.00] ; 
% Macro0{4}  = FACT*[-0.25; 0.10; 0.00; 0.00; 0.00] ;
% Macro0{5}  = FACT*[-0.20; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{6}  = FACT*[-0.06;-0.10; 0.00; 0.00; 0.00] ;
% Macro0{7}  = FACT*[ 0.00;-0.15; 0.00; 0.00; 0.00] ;
% Macro0{8}  = FACT*[ 0.19;-0.10; 0.00; 0.00; 0.00] ;
% Macro0{9}  = FACT*[ 0.00; 0.00; 0.00; 0.30; 0.00] ;
% Macro0{10} = FACT*[ 0.14; 0.00; 0.00; 0.33; 0.00] ;
% Macro0{11} = 1.60*[ 0.04; 0.04; 0.00; 0.10; 0.00] ;
% Macro0{12} = 1.50*[ 0.00; 0.07; 0.00; 0.10; 0.00] ;
% Macro0{13} = 3.80*[-0.05; 0.10; 0.00; 0.10; 0.00] ;
% Macro0{14} = 4.50*[-0.05; 0.00; 0.00; 0.10; 0.00] ;
% Macro0{15} = FACT*[-0.10;-0.10; 0.00; 0.15; 0.00] ;
% Macro0{16} = 2.50*[ 0.00;-0.10; 0.00; 0.20; 0.00] ;
% Macro0{17} = FACT*[ 0.19;-0.10; 0.00; 0.35; 0.00] ;

%% RVE_GrapInclu01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macro0{1}  = FACT*[ 0.10; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{2}  = FACT*[ 0.075; 0.075; 0.00; 0.00; 0.00] ; 
% Macro0{3}  = FACT*[ 0.00; 0.10; 0.00; 0.00; 0.00] ; 
% Macro0{4}  = FACT*[-0.25; 0.10; 0.00; 0.00; 0.00] ;
% Macro0{5}  = FACT*[-0.30; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{6}  = FACT*[-0.06;-0.10; 0.00; 0.00; 0.00] ;
% Macro0{7}  = FACT*[ 0.00;-0.15; 0.00; 0.00; 0.00] ;
% Macro0{8}  = FACT*[ 0.19;-0.10; 0.00; 0.00; 0.00] ;
% Macro0{9}  = FACT*[ 0.00; 0.00; 0.00; 0.30; 0.00] ;
% Macro0{10} = 0.97*[ 0.16; 0.00; 0.00; 0.37; 0.00] ;
% Macro0{11} = 1.60*[ 0.04; 0.04; 0.00; 0.10; 0.00] ;
% Macro0{12} = 1.50*[ 0.00; 0.07; 0.00; 0.10; 0.00] ;
% Macro0{13} = 3.80*[-0.05; 0.10; 0.00; 0.10; 0.00] ;
% Macro0{14} = 4.50*[-0.05; 0.00; 0.00; 0.10; 0.00] ;
% Macro0{15} = 1.90*[-0.07;-0.07; 0.00; 0.10; 0.00] ;
% Macro0{16} = 2.50*[ 0.00;-0.10; 0.00; 0.20; 0.00] ;
% Macro0{17} = 4.00*[ 0.05;-0.03; 0.00; 0.10; 0.00] ;

%% RVE_Hole01 Elas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macro0{1}  = FACT*[ 0.40; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{2}  = FACT*[ 0.30; 0.30; 0.00; 0.00; 0.00] ;
% Macro0{3}  = FACT*[ 0.00; 0.50; 0.00; 0.00; 0.00] ;
% Macro0{4}  = FACT*[-0.30; 0.40; 0.00; 0.00; 0.00] ;
% Macro0{5}  = FACT*[-0.25; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{6}  = FACT*[-0.15;-0.20; 0.00; 0.00; 0.00] ;
% Macro0{7}  = FACT*[ 0.00;-0.22; 0.00; 0.00; 0.00] ;
% Macro0{8}  = FACT*[ 0.40;-0.20; 0.00; 0.00; 0.00] ;
% Macro0{9}  = FACT*[ 0.00; 0.00; 0.00; 0.30; 0.00] ;
% Macro0{10} = FACT*[ 0.40; 0.00; 0.00; 0.30; 0.00] ;
% Macro0{11} = FACT*[ 0.30; 0.30; 0.00; 0.30; 0.00] ;
% Macro0{12} = FACT*[ 0.00; 0.40; 0.00; 0.30; 0.00] ;
% Macro0{13} = FACT*[-0.20; 0.50; 0.00; 0.30; 0.00] ;
% Macro0{14} = FACT*[-0.10; 0.00; 0.00; 0.30; 0.00] ;
% Macro0{15} = FACT*[-0.05;-0.05; 0.00; 0.30; 0.00] ;
% Macro0{16} = FACT*[ 0.00;-0.10; 0.00; 0.30; 0.00] ;
% Macro0{17} = FACT*[ 0.40;-0.20; 0.00; 0.30; 0.00] ;

%% RVE_Hole01 Plas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macro0{1}  = FACT*[ 0.14; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{2}  = FACT*[ 0.07; 0.07; 0.00; 0.00; 0.00] ;
% Macro0{3}  = FACT*[ 0.00; 0.08; 0.00; 0.00; 0.00] ;
% Macro0{4}  = FACT*[-0.08; 0.10; 0.00; 0.00; 0.00] ;
% Macro0{5}  = FACT*[-0.18; 0.00; 0.00; 0.00; 0.00] ;
% Macro0{6}  = FACT*[-0.05;-0.10; 0.00; 0.00; 0.00] ;
% Macro0{7}  = FACT*[ 0.00;-0.12; 0.00; 0.00; 0.00] ;
% Macro0{8}  = FACT*[ 0.20;-0.10; 0.00; 0.00; 0.00] ;
% Macro0{9}  = FACT*[ 0.00; 0.00; 0.00; 0.27; 0.00] ;
% Macro0{10} = FACT*[ 0.05; 0.00; 0.00; 0.10; 0.00] ;
% Macro0{11} = FACT*[ 0.04; 0.04; 0.00; 0.10; 0.00] ;
% Macro0{12} = FACT*[ 0.00; 0.07; 0.00; 0.10; 0.00] ;
% Macro0{13} = FACT*[-0.05; 0.10; 0.00; 0.10; 0.00] ;
% Macro0{14} = FACT*[-0.05; 0.00; 0.00; 0.10; 0.00] ;
% Macro0{15} = FACT*[-0.07;-0.07; 0.00; 0.10; 0.00] ;
% Macro0{16} = FACT*[ 0.00;-0.10; 0.00; 0.20; 0.00] ;
% Macro0{17} = FACT*[ 0.05;-0.03; 0.00; 0.10; 0.00] ;

%% RVE_Hole05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macro0{1}  = FACT*[0.07; 0; 0; 0; 0]; % Modo01  
% Macro0{2}  = FACT*[0; 0.07; 0; 0; 0]; % Modo02 
% Macro0{3}  = FACT*[0; 0; 0; 0.25; 0]; % Modo03 
% Macro0{4}  = FACT*[0.05; 0.05; 0; 0; 0]; % Modo04 
% Macro0{5}  = FACT*[0.07; 0; 0; 0.1; 0]; % Modo05
% Macro0{6}  = FACT*[-0.2; 0; 0; 0; 0]; % Modo06 
% Macro0{7}  = FACT*[-0.1; 0.05; 0; 0; 0.0]; % Modo07
% Macro0{8}  = FACT*[-0.1; 0; 0; 0.2; 0.0]; % Modo08
% Macro0{9}  = FACT*[-0.1; 0.1; 0; 0.3; 0]; % Modo09
% Macro0{10} = FACT*[-0.1; 0.1; 0; -0.2; 0]; % Modo10
% Macro0{11} = FACT*[0; 0.05; 0; -0.1; 0]; % Modo11
% Macro0{12} = FACT*[0; 0; 0; -0.3; 0]; % Modo12
% Macro0{13} = FACT*[0.05; 0.05; 0; -0.3; 0]; % Modo13
% Macro0{14} = FACT*[0.05; 0.05; 0; 0.1; 0]; % Modo14
% Macro0{15} = FACT*[0.04; 0; 0; -0.2; 0]; % Modo15
% Macro0{16} = FACT*[-0.05; 0; 0; -0.3; 0]; % Modo16 
% Macro0{17} = FACT*[0.07; 0; 0; 0.1; 0]; % Modo17
% Macro0{18} = FACT*[0; -0.12; 0; 0; 0]; % Modo18 
% Macro0{19} = FACT*[0; -0.12; 0; 0.2; 0]; % Modo19
% Macro0{20} = FACT*[0; -0.12; 0; -0.2; 0]; % Modo20
% Macro0{21} = FACT*[0.15; -0.15; 0; 0; 0]; % Modo21
% Macro0{22} = FACT*[0.15; -0.15; 0; 0.2; 0]; % Modo22
% Macro0{23} = FACT*[0.12; -0.12; 0; -0.2; 0]; % Modo23
% Macro0{24} = FACT*[-0.12; -0.12; 0; 0; 0]; % Modo24
% Macro0{25} = FACT*[-0.1; -0.1; 0; 0.2; 0]; % Modo25
% Macro0{26} = FACT*[0.1; -0.1; 0; -0.2; 0]; % Modo26

%% RVE_Struct_Hole_10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isMICRO(1).epsilon_Macro0=FACT*[0.2; 0; 0; 0; 0] ; %Modo01 ne11
% isMICRO(1).epsilon_Macro0=FACT*[-0.2; 0; 0; 0; 0] ; %Modo02 ne10
% isMICRO(1).epsilon_Macro0=FACT*[ 0; 0.2; 0; 0; 0] ; %Modo03 ne11
% isMICRO(1).epsilon_Macro0=FACT*[ 0; -0.2; 0; 0; 0] ; %Modo04 ne10
% isMICRO(1).epsilon_Macro0=FACT*[ 0; 0; 0; 0.2; 0.2] ; %Modo05 ne8
% isMICRO(1).epsilon_Macro0=FACT*[ 0; 0; 0;-0.2;-0.2] ; %Modo06 ne8
% isMICRO(1).epsilon_Macro0=FACT*[0.2; 0.2; 0; 0; 0] ; %Modo07 ne8
% isMICRO(1).epsilon_Macro0=FACT*[-0.14;-0.14; 0; 0; 0] ; %Modo08 ne17
% isMICRO(1).epsilon_Macro0=FACT*[0.2;-0.2; 0; 0; 0] ; %Modo09 ne8
% isMICRO(1).epsilon_Macro0=FACT*[-0.2; 0.2; 0; 0; 0] ; %Modo10 ne8
% isMICRO(1).epsilon_Macro0=FACT*[0.20; 0.0; 0.0; 0.2; 0.2] ; %Modo11 ne6
% isMICRO(1).epsilon_Macro0=FACT*[-0.10; 0.0; 0.0; 0.2; 0.2] ; %Modo12 ne8
% isMICRO(1).epsilon_Macro0=FACT*[0.0; 0.2; 0.0; 0.2; 0.2] ; %Modo13 ne6
% isMICRO(1).epsilon_Macro0=FACT*[0.0;-0.15; 0; 0.2; 0.2] ; %Modo14 ne7
% isMICRO(1).epsilon_Macro0=FACT*[0.2;-0.2; 0; 0.2; 0.2] ; %Modo15 ne6
% isMICRO(1).epsilon_Macro0=FACT*[-0.1;-0.1; 0; 0.2; 0.2] ; %Modo16 ne7
% isMICRO(1).epsilon_Macro0=FACT*[0.15; 0.15; 0; 0.2; 0.2] ; %Modo17 ne6
% isMICRO(1).epsilon_Macro0=FACT*[-0.2; 0.2; 0; 0.2; 0.2] ; %Modo18 ne6
% isMICRO(1).epsilon_Macro0=FACT*[0.2; 0.0; 0;-0.2;-0.2] ; %Modo19 ne6
% isMICRO(1).epsilon_Macro0=FACT*[-0.15; 0.0; 0;-0.2;-0.2] ; %Modo20 ne7
% isMICRO(1).epsilon_Macro0=FACT*[0.0; 0.2; 0;-0.2;-0.2] ; %Modo21 ne6
% isMICRO(1).epsilon_Macro0=FACT*[0.0;-0.15; 0;-0.2;-0.2] ; %Modo22 ne7
% isMICRO(1).epsilon_Macro0=FACT*[0.2;-0.2; 0;-0.2;-0.2] ; %Modo23 ne6
% isMICRO(1).epsilon_Macro0=FACT*[-0.1;-0.1; 0;-0.2;-0.2] ; %Modo24 ne6
% isMICRO(1).epsilon_Macro0=FACT*[0.15; 0.15; 0;-0.2;-0.2] ; %Modo25 ne6
% isMICRO(1).epsilon_Macro0=FACT*[-0.2; 0.2; 0;-0.2;-0.2] ; %Modo26 ne6
%% RVE_Rigid_Hole_Periodico %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isMICRO(1).epsilon_Macro0=FACT*[-0.23; 0; 0; 0; 0]; % Modo01   ne9 ne8 ne9
% isMICRO(1).epsilon_Macro0=FACT*[0; 0.15; 0; 0; 0]; % Modo02  ne14 ne13 ne9
% isMICRO(1).epsilon_Macro0=FACT*[0; 0; 0; 0.25; 0]; % Modo03  ne10 ne9 ne9 Hay que cambiar las funciones de carga en .mfl
% isMICRO(1).epsilon_Macro0=FACT*[0.15; -0.15; 0; 0; 0]; % Modo04  ne9 ne9 ne9
% isMICRO(1).epsilon_Macro0=FACT*[0.2; 0; 0; 0.1; 0]; % Modo05  ne10 ne9 ne10
% isMICRO(1).epsilon_Macro0=FACT*[0.1; -0.1; 0; 0.2; 0]; % Modo06 ne10 ne10 ne9
% isMICRO(1).epsilon_Macro0=FACT*[0.1; 0.1; 0; 0; -0.2]; % Modo07 ne9 ne8 ne9
% isMICRO(1).epsilon_Macro0=FACT*[-0.14; -0.14; 0; 0; 0]; % Modo08 ne12 ne11 ne9
% isMICRO(1).epsilon_Macro0=FACT*[-0.2; 0.1; 0; -0.1; -0.1]; % TRAYECTORIA NO ENTRENADA
%% RVE_Hole03 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isMICRO(1).epsilon_Macro0=FACT*[0.06; 0; 0; 0; 0]; % Modo01 ne22
% isMICRO(1).epsilon1_Macro0=FACT*[0; -0.06; 0; 0; 0]; % Modo02 ne23
% isMICRO(1).epsilon_Macro0=FACT*[0; 0; 0; 0.14; 0]; % Modo03 ne18
% isMICRO(1).epsilon_Macro0=FACT*[-0.06; 0.06; 0; 0; 0]; % Modo04 ne22
% isMICRO(1).epsilon_Macro0=FACT*[0.05; 0; 0; 0.05; 0]; % Modo05 ne22
% isMICRO(1).epsilon_Macro0=FACT*[0.05; -0.05; 0; 0.05; 0]; % Modo06 ne23
% isMICRO(1).epsilon_Macro0=FACT*[0.03; 0.03; 0; -0.05; 0.0]; % Modo07 ne24
% isMICRO(1).epsilon_Macro0=FACT*[-0.05; -0.05; 0; 0.0; 0.0]; % Modo08 ne25
% isMICRO(1).epsilon_Macro0=FACT*[0.07; 0; 0; 0.1; 0]; % Modo17
%% TEST 13 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isMICRO(1).epsilon_Macro0=FACT*[2.5; 0; 0; 0; 0]; % Modo01   ne10
% isMICRO(1).epsilon_Macro0=FACT*[0; -0.5; 0; 0; 0]; % Modo02    ne22
% isMICRO(1).epsilon_Macro0=FACT*[0; 0; 0; 2; 0]; % Modo03    ne8
% isMICRO(1).epsilon_Macro0=FACT*[0; 0; 0; 0; -2]; % Modo04    ne8
% isMICRO(1).epsilon_Macro0=FACT*[2.5; 2.5; 0; 0; 0]; % Modo05    ne10
% isMICRO(1).epsilon_Macro0=FACT*[1; 0; 0; 1; 0]; % Modo06  ne16
% isMICRO(1).epsilon_Macro0=FACT*[0; -0.7; 0; 0.5; 0]; % Modo07 ne14
% isMICRO(1).epsilon_Macro0=FACT*[1; -0.5; 0; 0.4; 0]; % Modo08 ne13
% isMICRO(1).epsilon_Macro0=FACT*[-0.5; 0.8; 0; 0; 0.3]; % TRAYECTORIA NO ENTRENADA
%% REV_Laminate_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isMICRO(1).epsilon_Macro0=FACT*[0.80; 0; 0; 0; 0]; % Modo01
% isMICRO(1).epsilon_Macro0=FACT*[0; 0.75; 0; 0; 0]; % Modo02 NO PLASTIFICA
% isMICRO(1).epsilon_Macro0=FACT*[0; 0; 0; 0; 1]; % Modo03 
% isMICRO(1).epsilon_Macro0=FACT*[0.7; -0.2; 0; 0; 0]; % Modo04
% isMICRO(1).epsilon_Macro0=FACT*[0.8; 0; 0; 1; 0]; % Modo05
% isMICRO(1).epsilon_Macro0=FACT*[0; 0.8; 0; 1; 0]; % Modo06
% isMICRO(1).epsilon_Macro0=FACT*[-0.5; 0.5; 0; 0; 0]; % Modo07
% isMICRO(1).epsilon_Macro0=FACT*[0.8; 0.8; 0; 1; 0]; % Modo08
% isMICRO(1).epsilon_Macro0=FACT*[0.8; -0.8; 0; 1; 1]; % Modo09
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
