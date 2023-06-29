#include "fintrf.h"

#if 0
C     
C     f_RMapPlastJ2LD.f
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif

C     [m_DPF,m_P,m_HVarNew,m_AuxVarOut] = f_RMapPlastJ2LD(m_F,m_HVarOld,m_AuxVarIn,e_DatMatSet,e_VG)
      subroutine mexFunction(nlhs,plhs,nrhs,prhs)


C     Se trat� de mantener compilable en Fortran 77.

C     Declarations
      implicit none
      
C     Constantes de dimensionamiento de arrays
      include 'ArraysDim.om'
      
C     Para suma de punteros de double
      integer pasoReal
      parameter(pasoReal=8)

C     mexFunction arguments:
      mwPointer plhs(*),prhs(*)
      integer nlhs,nrhs
      
C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxGetField
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateDoubleScalar
      real*8 mxGetScalar
      integer mxIsDouble
      integer mxIsStruct
      integer mexPrintf     
      
C     mwSize mxGetM, mxGetN
      
C     Pointers to input/output mxArrays:
      mwPointer p_F
      mwPointer p_HVarOld
      mwPointer p_HVarNew
C     mwPointer p_AuxVarIn
C     mwPointer p_AuxVarOut  
      mwPointer p_P
      mwPointer p_DPF

C     Variables auxiliares
      mwSize nTens       
      
C     Variables para transferencia y recepci�n de variables a la funci�n STRL81ld_matlab:
C      real*8 FVoigt[allocatable](:)
C     Deformation gradiente
      real*8 stran(stranDim,stranDim)
      real*8 stranzz         
C     Propiedades materiales      
      real*8 props(NPROP)
C     Tensor de Cauchy de deformaci�n derecho pl�stico (sim�trico)
C     CP = [Cxx,Cyy,Cxy,Czz]
      real*8 cplas(4)
C     Vector de tensiones de Kirchhoff (sim�trico) (tau = 
C     tau = [tauxx,tauyy,tauxy,tauzz]
      real*8 strsg(NSTR1)
C     Derivada derivada D[tau]/D[gradS[u(x)]]
C     Se asume la siguiente distribuci�n  [xx,yy,xy,zz]
      real*8 dmatx(NSTR1,NSTR1)
C     Deformaci�n pl�stica equivalente
      real*8 eqpst
C     Tensi�n equivalente pl�stica
      real*8 qfunr
C     Factor de carga
      real*8 fload
C     Variables para el paso de informaci�n entre la funci�n de c�lculo de tensiones y la funci�n del tensor
C     tangente constitutivo.
      real*8 factm,factp,factq
      real*8 vptria(vectorDim),vp(vectorDim)
      real*8 dirp1(vectorDim),dirp2(vectorDim)
      real*8 hpact
C     Matriz gradiente de deformaci�n para acceder desde fortran 
      real*8 F(nTensFijo)
C     Matriz de variables hist�ricas para acceder desde fortran
      real*8 HVar(nVarHist)
C     Determinante del gradiente de deformaci�n
      real*8 J
C     Matriz 2D inversa del gradiente de deformaci�n
      real*8 FInv2D(stranDim,stranDim)
C     T�rmino de zz de la inersa del gradiente de deformaci�n
      real*8 FInvzz
C     Primer tensor de tensiones de Piola-Kirchhoff traspuesto
      real*8 P(nTensFijo)
C     Tensor tangente constitutivo en notaci�n de Voigt (DP/DF)
C     Se asume esta distribuci�n [xx,yy,zz,xy]
      real*8 DPF(nTensFijo,nTensFijo)
C     Variables temporales para el c�lculo de DPF (para reducir el n�mero de operaciones)
      real*8 SqFInv11,SqFInv22,SqFInv33,SqFInv12,SqFInv21
C     Variables para calcular la energia
      real*8 Psi_e,Psi_ev,Psi_ed
      real*8 DTRES,YOUNG,POISS,MULAM,LAMBL2,KAPPA,trLogVP
      real*8 VPdev(vectorDim)
      
C     Variables para impresi�n 
C     integer ind,indi,indj
C     integer outPrint
C     character(len=100) textPrint
      
C     Check for proper number of arguments. 
      if(nrhs .lt. 5) then
         call mexErrMsgIdAndTxt('MEFNL:RMAPPlastJ2',
     1'Material Return Mapping: Plasticity 
     2J2 Large Deformation: Six input arguments is required.')
      elseif(nlhs .gt. 7) then
         call mexErrMsgIdAndTxt('MEFNL:RMAPPlastJ2',
     1'Material Return Mapping: Plasticity 
     2J2 Large Deformation: Too many output arguments.')
      endif
      
C     Validate inputs
C     Check that the first input is a number.
      if(mxIsDouble(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt('MEFNL:RMAPPlastJ2',
     1'Material Return Mapping: Plasticity 
     2J2 Large Deformation: First input must be a number matrix.')
      endif

C     Check that the second input is a number.      
      if(mxIsDouble(prhs(2)) .eq. 0) then
         call mexErrMsgIdAndTxt('MEFNL:RMAPPlastJ2',
     1'Material Return Mapping: Plasticity 
     2J2 Large Deformation: Second input must be a number matrix.')
      endif
      
C     Check that the third input is a number.      
      if(mxIsDouble(prhs(3)) .eq. 0) then
         call mexErrMsgIdAndTxt('MEFNL:RMAPPlastJ2',
     1'Material Return Mapping: Plasticity 
     2J2 Large Deformation: Third input must be a number matrix.')
      endif
      
C     Check that the fourth is a struct.
      if(mxIsStruct(prhs(4)) .eq. 0) then
         call mexErrMsgIdAndTxt('MEFNL:RMAPPlastJ2',
     1'Material Return Mapping: Plasticity 
     2J2 Large Deformation: Fourth input must be a struct.')
      endif
      
C     Check that the fifth is a struct.
      if(mxIsStruct(prhs(5)) .eq. 0) then
         call mexErrMsgIdAndTxt('MEFNL:RMAPPlastJ2',
     1'Material Return Mapping: Plasticity 
     2J2 Large Deformation: Fifth input must be a struct.')
      endif
      
C     Componentes de tensi�n y deformaci�n dentro del programa de matlab
C     nTens = mxGetScalar(mxGetField(prhs(5),1,'ntens'))
      
C     Impresi�n en l�nea de comando 
C      write(textPrint,*) 'nTens = ',nTens
C      outPrint = mexPrintf(textPrint//achar(10))
      
C     Gradient Deformation Tensor in fortran format para transferencia a la funci�n
      p_F = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(p_F,F,nTensFijo)
C     stran(1,1) = mxGetScalar(prhs(1))
      stran(1,1) = F(1)
      stran(2,2) = F(2)
      stran(1,2) = F(4)
      stran(2,1) = F(5)
C     Impresi�n en l�nea de comando
C     do indi = 1,stranDim
C        do indj = 1,stranDim
C           write(textPrint,*) 'stran ',stran(indi,indj)
C           outPrint = mexPrintf(textPrint//achar(10))
C        enddo
C     enddo
      stranzz = F(3)
C     Impresi�n en l�nea de comando 
C     write(textPrint,*) 'stranzz = ',stranzz
C     outPrint = mexPrintf(textPrint//achar(10))

C     Vector de propiedades materiales      
      props(1) = mxGetScalar(mxGetField(prhs(4),1,'young'))
      props(2) = mxGetScalar(mxGetField(prhs(4),1,'poiss'))
      props(3) = mxGetScalar(mxGetField(prhs(4),1,'yield'))
      props(4) = mxGetScalar(mxGetField(prhs(4),1,'hba'))
      props(5) = mxGetScalar(mxGetField(prhs(4),1,'tit'))
      props(6) = mxGetScalar(mxGetField(prhs(4),1,'kin'))
      props(7) = mxGetScalar(mxGetField(prhs(4),1,'kce'))
      props(8) = mxGetScalar(mxGetField(prhs(4),1,'del'))
C     Impresi�n en l�nea de comando
C     do ind = 1,NPROP
C        write(textPrint,*) 'props ',props(ind)
C        outPrint = mexPrintf(textPrint//achar(10))
C     enddo
      
C     Extraci�n de variables hist�ricas
      p_HVarOld = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(p_HVarOld,HVar,nVarHist)
C     Variable deformaci�n pl�stica equivalente (se asume en el primer valor de m_HVarOld)
      eqpst = HVar(1)
C     Variable tensi�n pl�stica equivalente (segundo valor)
      qfunr = HVar(2)
C     Tensor de Cauchy de deformaci�n derecho pl�stico (sim�trico)
C     Componentes XX
      cplas(1) = HVar(3)
C     Componente YY
      cplas(2) = HVar(4)
C     Componente XY 
      cplas(3) = HVar(6)
C     Componente ZZ
      cplas(4) = HVar(5)
      
C     Extracci�n de variables auxiliares
C     p_AuxVar = mxGetPr(prhs(3))
C     Factor de carga (como por ahora hay una sola se utiliza directamente mxGetScalar)
      fload = mxGetScalar(prhs(3))

C     Determinaci�n de las tensiones de Kirchoff tau (tau = J * Sigma, Sigma: Tensi�n de Cauchy)
      call STRL81ld_matlab(props,stran,strsg,eqpst,qfunr,fload,
     1factm,factp,factq,
     2hpact,
     3cplas,
     4vptria,vp,dirp1,dirp2,
     5stranzz)
     
C     Impresi�n en l�nea de comando
C     do ind = 1,NSTR1
C        write(textPrint,*) 'tens ',strsg(ind)
C        outPrint = mexPrintf(textPrint//achar(10))
C     enddo

C     Determinaci�n de la derivada D[tau]/D[gradS[u(x)]]  
      call MTDL81(stran,props,dmatx,fload,
     1factm,factp,factq,
     2vptria,vp,dirp1,dirp2)
      
C     Determinante del gradiente de deformaci�n (J = det(F))
      J = (stran(1,1)*stran(2,2)-stran(1,2)*stran(2,1))*stranzz

C     Inversa del gradiente de deformaci�n
      FInv2D(1,1) = F(2)/J
      FInv2D(2,2) = F(1)/J
      FInv2D(1,2) = -F(4)/J
      FInv2D(2,1) = -F(5)/J
      FInvzz = 1/F(3)      
      
C     C�lculo de Primer tensor de tensiones de Piola-Kirchhoff
C     Componente xx
      P(1) = strsg(1)*FInv2D(1,1)+strsg(3)*FInv2D(1,2)
C     Componente yy
      P(2) = strsg(3)*FInv2D(2,1)+strsg(2)*FInv2D(2,2)
C     Componente zz
      P(3) = strsg(4)*FInvzz
C     Componente xy 
      P(4) = strsg(1)*FInv2D(2,1)+strsg(3)*FInv2D(2,2)
C     Componente yx      
      P(5) = strsg(3)*FInv2D(1,1)+strsg(2)*FInv2D(1,2)
      
C     Como se utiliza en varias t�rminos se calcula los cuadrados de los t�rminos de FInv
      SqFInv11 = FInv2D(1,1)**2
      SqFInv22 = FInv2D(2,2)**2
      SqFInv33 = FInvzz**2
      SqFInv12 = FInv2D(1,2)**2
      SqFInv21 = FInv2D(2,1)**2
      
C     Determinaci�n de la derivada D[P]/D[F]
C     Componente 1111
      DPF(1,1) = SqFInv11*(dmatx(1,1)+strsg(1))+
     1SqFInv12*(dmatx(3,3)+strsg(2))+
     2FInv2D(1,1)*FInv2D(1,2)*(dmatx(1,3)+dmatx(3,1)+2*strsg(3))
C     Componente 2211
      DPF(2,1) = FInv2D(1,1)*(dmatx(3,1)*FInv2D(2,1)+
     1dmatx(2,1)*FInv2D(2,2))+
     2FInv2D(1,2)*(dmatx(3,3)*FInv2D(2,1)+dmatx(2,3)*FInv2D(2,2))
C     Componente 3311
      DPF(3,1) = (dmatx(4,1)*FInv2D(1,1)+dmatx(4,3)*FInv2D(1,2))*FInvzz
C     Componente 1211
      DPF(4,1) = FInv2D(1,2)*(FInv2D(2,2)*(dmatx(3,3)+strsg(2))+
     1FInv2D(2,1)*(dmatx(1,3)+strsg(3)))+
     2FInv2D(1,1)*(FInv2D(2,1)*(dmatx(1,1)+strsg(1))+
     3FInv2D(2,2)*(dmatx(3,1)+strsg(3)))
C     Componente 2111
      DPF(5,1) = (dmatx(2,1)+dmatx(3,3))*FInv2D(1,1)*FInv2D(1,2)+
     1dmatx(3,1)*SqFInv11+dmatx(2,3)*SqFInv12
C     Componente 1122
      DPF(1,2) = FInv2D(1,1)*(dmatx(1,3)*FInv2D(2,1)+
     1dmatx(1,2)*FInv2D(2,2))+
     2FInv2D(1,2)*(dmatx(3,3)*FInv2D(2,1)+
     3dmatx(3,2)*FInv2D(2,2))
C     Componente 2222
      DPF(2,2) = SqFInv21*(dmatx(3,3)+strsg(1))+
     1SqFInv22*(dmatx(2,2)+strsg(2))+
     2FInv2D(2,1)*FInv2D(2,2)*(dmatx(2,3)+dmatx(3,2)+2*strsg(3))
C     Componente 3322
      DPF(3,2) = (dmatx(4,3)*FInv2D(2,1)+dmatx(4,2)*FInv2D(2,2))*FInvzz
C     Componente 1222
      DPF(4,2) = (dmatx(1,2)+dmatx(3,3))*FInv2D(2,1)*FInv2D(2,2)+
     1dmatx(1,3)*SqFInv21+dmatx(3,2)*SqFInv22
C     Componente 2122     
      DPF(5,2) = FInv2D(1,2)*(FInv2D(2,2)*(dmatx(2,2)+strsg(2))+
     1FInv2D(2,1)*(dmatx(2,3)+strsg(3)))+
     2FInv2D(1,1)*(FInv2D(2,1)*(dmatx(3,3)+strsg(1))+
     3FInv2D(2,2)*(dmatx(3,2)+strsg(3)))
C     Componente 1133
      DPF(1,3) = (dmatx(1,4)*FInv2D(1,1)+dmatx(3,4)*FInv2D(1,2))*FInvzz
C     Componente 2233
      DPF(2,3) = (dmatx(3,4)*FInv2D(2,1)+dmatx(2,4)*FInv2D(2,2))*FInvzz
C     Componente 3333
      DPF(3,3) = SqFInv33*(dmatx(4,4)+strsg(4))
C     Componente 1233
      DPF(4,3) = (dmatx(1,4)*FInv2D(2,1)+dmatx(3,4)*FInv2D(2,2))*FInvzz
C     Componente 2133
      DPF(5,3) = (dmatx(3,4)*FInv2D(1,1)+dmatx(2,4)*FInv2D(1,2))*FInvzz
C     Componente 1112
      DPF(1,4) = FInv2D(1,1)*(FInv2D(2,1)*(dmatx(1,1)+strsg(1))+
     1FInv2D(2,2)*(dmatx(1,3)+strsg(3)))+
     2FInv2D(1,2)*(FInv2D(2,2)*(dmatx(3,3)+strsg(2))+
     3FInv2D(2,1)*(dmatx(3,1)+strsg(3)))
C     Componente 2212
      DPF(2,4) = (dmatx(2,1)+dmatx(3,3))*FInv2D(2,1)*FInv2D(2,2)+
     1dmatx(3,1)*SqFInv21+
     2dmatx(2,3)*SqFInv22
C     Componente 3312
      DPF(3,4) = (dmatx(4,1)*FInv2D(2,1)+dmatx(4,3)*FInv2D(2,2))*FInvzz
C     Componente 1212
      DPF(4,4) = SqFInv21*(dmatx(1,1)+strsg(1))+
     1SqFInv22*(dmatx(3,3)+strsg(2))+
     2FInv2D(2,1)*FInv2D(2,2)*(dmatx(1,3)+dmatx(3,1)+2*strsg(3))
C     Componente 2112
      DPF(5,4) = FInv2D(1,2)*(dmatx(2,1)*FInv2D(2,1)+
     1dmatx(2,3)*FInv2D(2,2))+
     2FInv2D(1,1)*(dmatx(3,1)*FInv2D(2,1)+dmatx(3,3)*FInv2D(2,2))
C     Componente 1121
      DPF(1,5) = (dmatx(1,2)+dmatx(3,3))*FInv2D(1,1)*FInv2D(1,2)+
     1dmatx(1,3)*SqFInv11+dmatx(3,2)*SqFInv12
C     Componente 2221
      DPF(2,5) = FInv2D(1,1)*(FInv2D(2,1)*(dmatx(3,3)+
     1strsg(1))+FInv2D(2,2)*(dmatx(2,3)+strsg(3)))+
     2FInv2D(1,2)*(FInv2D(2,2)*(dmatx(2,2)+strsg(2))+
     3FInv2D(2,1)*(dmatx(3,2)+strsg(3)))
C     Componente 3321
      DPF(3,5) = (dmatx(4,3)*FInv2D(1,1)+dmatx(4,2)*FInv2D(1,2))*FInvzz
C     Componente 1221
      DPF(4,5) = FInv2D(1,2)*(dmatx(1,2)*FInv2D(2,1)+
     1dmatx(3,2)*FInv2D(2,2))+
     2FInv2D(1,1)*(dmatx(1,3)*FInv2D(2,1)+
     3dmatx(3,3)*FInv2D(2,2))
C     Componente 2121
      DPF(5,5) = SqFInv11*(dmatx(3,3)+strsg(1))+
     1SqFInv12*(dmatx(2,2)+strsg(2))+
     2FInv2D(1,1)*FInv2D(1,2)*(dmatx(2,3)+dmatx(3,2)+2*strsg(3))
      
C     Se crea la matrix mxMatrix para devolver como argumento el tensor tangente constitutivo
      plhs(1) = mxCreateDoubleMatrix(nTensFijo,nTensFijo,0)
      p_DPF = mxGetPr(plhs(1))
C     Se copia a la matriz de MatLab
      call mxCopyReal8ToPtr(DPF,p_DPF,nTensFijo*nTensFijo)
      
C     Se crea la matrix mxMatrix para devolver como argumento la tensi�n
      plhs(2) = mxCreateDoubleMatrix(nTensFijo,1,0)
      p_P = mxGetPr(plhs(2))
C     Se copia a la matriz de MatLab
      call mxCopyReal8ToPtr(P,p_P,nTensFijo)
      
C     Se obtiene una array de fortran ordenado seg�n la convenci�n usada en el
C     programa de MatLab para la variable hist�rica new.
      HVar(1) = eqpst
      HVar(2) = qfunr
      HVar(3) = cplas(1)
      HVar(4) = cplas(2)
      HVar(5) = cplas(4) 
      HVar(6) = cplas(3)


      DTRES = 2.D0/3.D0
      YOUNG = PROPS(1)                      ! YOUNG MODULUS
      POISS = PROPS(2)                      ! POISSON RATIO
C     MULAN: SHEAR MODULUS (G)
      MULAM = YOUNG/2.D0/(1.D0+POISS)
C     LAMBL: LAME'S FIRST PARAMETER
      LAMBL2=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)/2.D0
C     KAPPA: BULK MODULUS (K) 2*LAMBL2+DTRES*MULAM
      KAPPA = YOUNG/3.D0/(1.D0-2.D0*POISS)

c     Elastic strain energy
      Psi_e = MULAM*((DLOG(VP(1)))**2 + (DLOG(VP(2)))**2 + 
     .      (DLOG(VP(3)))**2 ) + LAMBL2*(DLOG(VP(1)*VP(2)*VP(3)))**2 
      trLogVP  = DLOG(VP(1)) + DLOG(VP(2)) + DLOG(VP(3))
      VPdev(1) = 2.D0/3.D0*DLOG(VP(1)) 
     .                        - 1.D0/3.D0*(DLOG(VP(2))+DLOG(VP(3)))
      VPdev(2) = 2.D0/3.D0*DLOG(VP(2)) 
     .                        - 1.D0/3.D0*(DLOG(VP(1))+DLOG(VP(3)))
      VPdev(3) = 2.D0/3.D0*DLOG(VP(3)) 
     .                        - 1.D0/3.D0*(DLOG(VP(2))+DLOG(VP(1)))

      Psi_ev = 1.D0/2.D0*KAPPA*trLogVP**2
      Psi_ed = MULAM*(VPdev(1)**2 + VPdev(2)**2 + VPdev(3)**2)
      HVar(7) = Psi_e 
      HVar(8) = Psi_ev
      HVar(9) = Psi_ed

C     Se crea una nueva variable mxMatrix para pasar la variable hist�rica new
C     (por si acaso no se usa el p_HVarOld, ya que capaz modifica la variable original,
C     ver como funciona eso)
      plhs(3) = mxCreateDoubleMatrix(nVarHist,1,0)
      p_HVarNew = mxGetPr(plhs(3))
C     Se copia a la matriz de MatLab
      call mxCopyReal8ToPtr(HVar,p_HVarNew,nVarHist)
      
C     Se crea una nueva variable mxMatrix para pasar la variable auxiliar
C     (Habr�a que ver si se puede usar la misma de entrada)
C     Como se tiene una sola variable auxiliar se usar la funci�n mxCreateDoubleScalar
      plhs(4) = mxCreateDoubleScalar(fload)
      
C     Se env�a para debug el tensor de tensiones de Kirchhoff (tau) y el tensor tangente 
C     D[tau]/D[gradS[u(x)]].
C     Tau:
      plhs(5) = mxCreateDoubleMatrix(NSTR1,1,0)
C     Se copia a la matriz de MatLab
      call mxCopyReal8ToPtr(strsg,mxGetPr(plhs(5)),NSTR1)
C     D[tau]/D[gradS[u(x)]]:
      plhs(6) = mxCreateDoubleMatrix(NSTR1,NSTR1,0)
C     Se copia a la matriz de MatLab
      call mxCopyReal8ToPtr(dmatx,mxGetPr(plhs(6)),NSTR1*NSTR1)

      return
      
      end

