% *************************************************************************
% MULTISCALE PROGRAM - MAIN FILE
% *************************************************************************
clear all; close all
main_example_path = pwd ;%'C:\Users\Hp 15-K118\Desktop\Projects\Multiscale_Application';%/home/javiermro/Projects/Examples'; 
%% ********* TESTS **********
first_mode=1;   nmode = 1;   Macro0 = cell(nmode,1);  
%#######################################################################################
%Para analisis MULTIESCALA
TEST_DATA(1).path_file= [main_example_path '/Ejemplos/RVE_1_Tay_Per_15_Analitico/']; TEST_DATA(1).file = 'Macro.mfl' ;
isMICRO.MICRO =0; % For MS analysis

%#######################################################################################
%Para analisis RVE
% TEST_DATA(1).path_file= [main_example_path '/Ejemplos/Microescala_RVE_64_elem/Test_RVE_17/']; TEST_DATA(1).file = 'RVE_Inclu01.mfl' ;
% isMICRO.MICRO =1; % For RVE analysis

%Ejemplo 1
% Elemento 5. PG 1. Step 1
% FACT = 1; Macro0{1}  = FACT*[-2.083505642557217e 1.641835190188765e-07 5.162263291583062e-12 0.0]'; 
% FACTp = 1; poreMacro0{1}  = FACTp*0.100030237480684'; 
% FACTphi = 1; gradPMacro0{1}  = FACTphi*[-8.005957280626096e-17 -0.002469982247133]'; 

%Ejemplo 2
% Elemento 5. PG 1. Step 12
% FACT = 1; Macro0{1}  = FACT*[-1.559357096203169e-18 -2.561786962506063e-06 5.162262776075984e-11 0.0]'; 
% FACTp = 1; poreMacro0{1}  = FACTp*0.999736455434032'; 
% FACTphi = 1; gradPMacro0{1}  = FACTphi*[-2.594596100937441e-15 0.0102014455968017]'; 

%Ejemplo 3
% Elemento 5. PG 1. Step 51
% FACT = 1; Macro0{1}  = FACT*[1.287663183787759e-18 -0.007428496354070 5.162262947483576e-11 0.0]'; 
% FACTp = 1; poreMacro0{1}  = FACTp*8.106858058018100e-06'; 
% FACTphi = 1; gradPMacro0{1}  = FACTphi*[-8.850427870901217e-18 -5.327702453714648e-06]'; 
%#######################################################################################

TEST_DATA(1).nLab = 0; % Parallel pools

% first_mode=9;   nmode = 9;   Macro0 = cell(nmode,1);
% for imodo=first_mode:nmode
%     if imodo<10
%         TEST_DATA(imodo).path_file= [main_example_path '/RVE_MetallicAlloy3/Modo0' num2str(imodo)]; 
%     else
%         TEST_DATA(imodo).path_file= [main_example_path '/RVE_MetallicAlloy3/Modo' num2str(imodo)]; 
%     end 
%     TEST_DATA(imodo).file = 'RVE_MetallicAlloy4000.mfl' ; 
%     TEST_DATA(imodo).nLab = 1;   
% end
% isMICRO.MICRO =1; FACT = 1.00;
% Trajectories ;

%% Snapshots que se desean calcular
SnapStrain       = false;
SnapStress       = false;
SnapWeight       = false;
SnapEnergy_e     = false;
SnapEnergy_e_vol = false;
SnapEnergy_e_dev = false;
SnapEnergy_p     = false;
SnapEnergy_t     = false;
Snapflag         = false;
% Snaps=[SnapStrain SnapStress SnapWeight SnapEnergy_e SnapEnergy_e_vol ...
%     SnapEnergy_e_dev SnapEnergy_p SnapEnergy_t Snapflag];
Snaps=[];
%%

% for iTEST = 1:length(TEST_DATA)
for iTEST=first_mode:nmode
    clc
    disp(['*** TEST NUMBER: ' num2str(iTEST) ' ***']);
%   try   
    isMICRO.epsilon_Macro0 = Macro0{iTEST};
%     isMICRO.p_Macro0 = poreMacro0{iTEST};
%     isMICRO.phi_Macro0 = gradPMacro0{iTEST};
    
    analysis(TEST_DATA(iTEST).path_file,TEST_DATA(iTEST).file,TEST_DATA(iTEST).nLab,isMICRO,Snaps);
    
%        matlabpool close
%   catch
%       warning('El proceso ha parado su ejecucion por falta de convergencia!');
%       matlabpool close
%   end
%
end

