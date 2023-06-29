function u = f_porp_nint(u,e_DatSet,e_VG)
% Permite determinar por interpolacion el valor de las

nSet = e_VG.nSet;
ndn = e_VG.ndn;

% Inicializaciones
% c_U = cell(nSet,1);
% c_FilU = cell(nSet,1);

   for iSet = 1:nSet
       nElem = e_DatSet(iSet).nElem;
       m_DofElem = e_DatSet(iSet).m_DofElem;
       conec = e_DatSet(iSet).conec;
       
       dofElemSet = m_DofElem(:);
%        m_FilU = dofElemSet';
       dofpe = e_DatSet(iSet).e_DatElem.dofpe;
       npe = e_DatSet(iSet).e_DatElem.npe;
       ndn_p = e_DatSet(iSet).e_DatElem.ndn_p;
       dofElemSet = m_DofElem(:);
       uElemSet  = reshape(u(dofElemSet),[],nElem);
       
       conec_int=conec(:,4*ndn_p+1:npe);
       m_FF_p = sf_nodmid(ndn_p,npe);
       
       pore_press = uElemSet(3:ndn:dofpe,:);
       con_ant = [];
       for iElem=1:nElem
           con_rep=[];
           con=conec_int(iElem,:);
           if iElem>1
               for icon=1:iElem-1
                   a1=repmat(con,4,1)';
                   a2=repmat(con_ant(icon,:),4,1);
                   con_id=(a1==a2);
                   a3=int8(con_id);
                   a4=sum(a3,2)';
                   for i=1:4
                       if a4(i)==1
                           con_rep=[con_rep 4*ndn_p+i];
                       end
                   end
               end
           end
           pp=pore_press(:,iElem);
           pore_press(4*ndn_p+1:npe,iElem)=m_FF_p*pp;
           pore_press(con_rep,iElem)=0.0;
           u(3*con)= u(3*con) + pore_press(4*ndn_p+1:npe,iElem);
%            u(3*con_rep)=0.0;
           con_ant=[con_ant ; con];
%            uElemSet(3:ndn:dofpe,iElem)=pore_press(1:npe,iElem);
       end %iElem
%        c_FilU{iSet} = m_FilU;
%        c_U{iSet} = uElemSet(:);
   end
%    % Ensamble de matriz de desplazamientos con poropresiones en lados del
%    % elemento
%    u = sparse([c_FilU{:}],1,cat(1,c_U{:}));
%    u=full(u);
     
   function m_FF_p = sf_nodmid(ndn_p,npe)
%*  numeraci�n de nodos: 4------7------3                                         *
%*                       |             |                                         *
%*                       8             6                                         *
%*                       |             |                                         *
%*                       1------5------2                                         *
   
   m_FF_p = zeros(4*ndn_p,npe); % Poropresiones 
   
   xg=[0 -1;1 0;0 1;-1 0]; %Cordenadas del centro de cada lado (Nodos 5, 6, 7 y 8)
   
   for ip=1:4
       E = xg(ip,1);
       n = xg(ip,2);
       % Funciones de forma: Elemento cuadrilatero de 4 nodos (P/poropresiones)
       SF1 = 1/4*(1-E).*(1-n);
       SF2 = 1/4*(1+E).*(1-n);
       SF3 = 1/4*(1+E).*(1+n);
       SF4 = 1/4*(1-E).*(1+n);
   
       % P/Poropresiones
       m_FF_p(ip,1)   = SF1; 
       m_FF_p(ip,2)   = SF2;
       m_FF_p(ip,3)   = SF3;
       m_FF_p(ip,4)   = SF4;
   end