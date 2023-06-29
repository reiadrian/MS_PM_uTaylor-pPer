clc
clear all
% load('SAVE.mat','y')
q = ones(10,1);
% for i=1:10    
    
    save('SAVE.mat','q','-append') ;
    q2=q*2 ;
    save('SAVE.mat','q2','-append') ;
% end
q
clear all
load('SAVE.mat','q')
q

% delete('/mytests/*.mat' %borra todos los archivos con extension .mat



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verificacion de Snapshots: integral de SnapStrain*DetJ*Wgauss = 0
% Verificacion de Modos: integral de Phi*DetJ*Wgauss = 0

4 * 513
W=[] ; for i=1:513 ;     W=[W ; e_DatSet.m_DetJT(:,i)]; end

W=ones(25*4,1)*0.01 ;


sum = 0;
ii = 1 ;
for ig=1:25*4
    for it=1:4    
        sum = sum + SnapStrain(ii,1)*W(ig) ;
        ii = ii+1 ;
    end
end
sum

sum = 0;
ii = 1 ;
for ig=1:513*4
    for it=1:4    
        sum = sum + PHI_GEN(ii,2)*W(ig) ;
        ii = ii+1 ;
    end
end