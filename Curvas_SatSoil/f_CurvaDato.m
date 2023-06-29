function tidespor = f_CurvaDato
% dp = f_CurvaDato(fid)

%*******************************************************************************
%* APERTURA DEL ARCHIVO DE DATOS                                               *
%*******************************************************************************
main_example_path =[pwd '/Examples/Articulo2/MS_RVE8/'];
path_file= main_example_path; file = 'CurvaDato.dat' ; %'\Examples\MS_Comp\'
file_open = fullfile(path_file,file);
fid = fopen(file_open,'r');

% Lectura de datos de tiempo, desplazamientos y poropresiones
nSteps = 51;
format = '%q %q %q';
tidespor = textscan(fid,format,nSteps,'CollectOutput',1,'CommentStyle','$');
tidespor = tidespor{1};
tidespor = str2double(tidespor(:,:,:));

%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);
