function [du_iter,dlambda] = incremental_disp(res,KT,Fext,m_LinCond,dofl,doff,Du_step_new,...
   Du_step_old,e_VG)

%******************************************************************************************
%*  COMPUTO DEL INCREMENTO DE DESPLAZAMIENTO ITERATIVO                                    *
%*  SEGUN LA ESTRATEGIA DE CONTROL SELECCIONADA                                           *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

CONTROL_STRAT = e_VG.CONTROL_STRAT;

% ********************************
% * INCREMENTO DE DESPLAZAMIENTO *
% ********************************
switch CONTROL_STRAT
   case 1
      KT = KT(dofl,dofl)+m_LinCond'*KT(doff,dofl)+KT(dofl,doff)*m_LinCond+...
         m_LinCond'*KT(doff,doff)*m_LinCond;
      %Considerando que la matriz es simetrica (se la fuerza para que rdivide la vea como simetrica
      %aun con los errores de precision).
      if e_VG.esKTSim
         KT = (KT+KT')/2;
      end
      du_iter = -KT\res;
      dlambda = 0.0;
   case 3
      [du_iter,dlambda,iconv] = update_normal_plane(res,KT,Fext,dofl,Du_step_new,Du_step_old,e_VG);
   case 4       
      Fext_uimp = KT(dofl,doff)* e_VG.vfix(doff) - Fext(dofl);
      % 
      KT = KT(dofl,dofl)+m_LinCond'*KT(doff,dofl)+KT(dofl,doff)*m_LinCond+...
         m_LinCond'*KT(doff,doff)*m_LinCond;
      %Considerando que la matriz es simetrica (se la fuerza para que rdivide la vea como simetrica
      %aun con los errores de precision).
      if e_VG.esKTSim
         KT = (KT+KT')/2;
      end
      [du_iter,dlambda,iconv] = control_desplazamientos(res,KT,Fext_uimp,e_VG.vfix,...
          dofl,doff,Du_step_new,Du_step_old,e_VG);
   otherwise
      error('Newton-Raphson: Incremento de desplazamiento: Estrategia de control incorrecta.')
end