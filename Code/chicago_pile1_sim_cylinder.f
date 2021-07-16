! *** Monte Carlo Simulation for Simplified Chicago Pile 1
!     Find Optimal Distance for Moderation of Fission Neutrons
!     In Graphite Plate
!     Neglect Effect of V of Target Nuclei
!     Epithermal: Neutron E<500 EV

      Program MC Pile1  
      Common IX,ANUMBE
      Dimension ENEOUT(100)
      IX=77777
      Open(1,File='Pile1.out',Status='UNKNOWN')
      Open(2,File='eneout.out',Status='UNKNOWN')


! *** Track BACKSCATTERED, TRANSMITTED, ABSORBED, THERMAL NEUTRONS
      BCKSCT=0
      TRANSM=0
      ABSORB=0
      THERMA=0

! *** Initialize Spectrum of Transmitted Energies (100 Bins)
      DO 3 J=1,100
      ENEOUT(J)=0
    3 CONTINUE

! *** Set Up Graphite Slab of Thickness D
      DENSIT=2.1
      ANUMBE=12.
      AVOGAD=6.02E23
      D=20.0
   50 Print*, "Enter Graphite Thickness D (in cm):"
      Read(*,*)D
      D = D
      
      BCKSCT=0
      TRANSM=0
      ABSORB=0
      THERMA=0

! *** Initialize Origin of Neutrons
      !Starting Histories of Neutrons
      RMC=1.E4
      DO 99 RX=1,RMC
      X=0
      Y=0
      Z=0
      
! *** Call Initial Energy Spectrum
      Call ENERGY(ENE)
      
! *** Randomize Entry Angle
      Call RANDM(RN1)
      ANGL=ACOS(RN1)
      Call RANDM(RN1)
      AZIM=2.*3.14159*RN1
      
! *** Build Velocity Unit Vectors
      VX=SIN(ANGL)*COS(AZIM)
      VY=SIN(ANGL)*SIN(AZIM)
      VZ=COS(ANGL)
      
! *** Call Cross Section Subroutine to Det Distance to Next Interaction -> DNEXT
!     Use Inverse Trasform of Probability Distribution l=-(ln(rnd#)*A)/(sigma*roe*Na)
    1 Call XSECT(ENE,XELAST,XABSOR)
      XTOTAL=XELAST+XABSOR
      Call RANDM(RN1)
      DNEXT=(-LOG(RN1)*ANUMBE)/(XTOTAL*DENSIT*AVOGAD)

! *** Calc New Position Using DNEXT (distance to next interaction)
      X=X+(VX*DNEXT)
      Y=Y+(VY*DNEXT)
      Z=Z+(VZ*DNEXT)
      
! *** Tracking Time -> Ask the Right Questions

! *** Has it BackScattered?
      If(Z.LT.0)Then
      BCKSCT=BCKSCT+1
      GOTO 99
      End If
! *** Has it Transmitted and if so was it Thermalized?
      If(Z.GT.D)Then
      TRANSM=TRANSM+1
      INDEX=INT(ENE*10.)+1
      ENEOUT(INDEX)=ENEOUT(INDEX)+1.
      If(ENE.LT.5E-4)THERMA=THERMA+1.
      GOTO 99
      End If
! *** Has it been Absorbed?
      RATIO=XABSOR/XTOTAL
      Call RANDM(RN1)
      If(RN1.LT.RATIO)Then
      ABSORB=ABSORB+1
      GOTO 99
      End If      

! *** Call ELOSS for Scattered Neutrons still Remaining
      Call ELOSS(ENE,ANGL,ENEAFT)
      ENE=ENEAFT
      
! *** Call EULER angles to find the new trajectory after ELOSS collision
      Call EULER(VX,VY,VZ,ANGL,SX,SY,SZ)
      VX=SX
      VY=SY
      VZ=SZ
      GOTO 1

   99 Continue
      
      !Percent=100./RMC
      Print*,RX
      Write(*,*)"Thickness D:",D,"cm"
      Write(*,*)"%TRANSM. %BCKSCT. %ABSORB. %THERMA."
      Write(*,*)TRANSM/RMC,BCKSCT/RMC,ABSORB/RMC,THERMA/RMC
      Write(1,*)"Thicknes D:",D,"cm"
      Write(1,*)"%TRANSM. %BCKSCT. %ABSORB. %THERMA."
      Write(1,*)TRANSM,BCKSCT,ABSORB,THERMA
      DO 20 G=1,100
      Write(2,*)G/10,ENEOUT(INT(G))
   20 Continue
      GOTO 50
      STOP 5555      
      End Program MC Pile1

! ************************************************
! ************************************************
      SUBROUTINE ELOSS(ENE,ANGL,ENEAFT)
      COMMON IX,ANUMBE
! *** CALCULATE ENERGY LOSS IN COLLISION, USE ISOTROPIC (S-WAVE) 
! *** APPROXIMATION. ROUND UP ANY ENERGY BELOW 0.1 EV TO 0.02 EV (THERMAL)
      IF(ENE.LT.1.E-7)THEN
      ENEAFT=2.E-8
      ELSE
! *** THIS IS THE MAXIMUM RECOIL ENERGY THAT CAN BE IMPARTED
      EREMAX=((4.*ANUMBE)/((1.+ANUMBE)**2.))*ENE
      CALL RANDM(RN1)
! *** THE PROBABILITY DISTRIBUTION FOR RECOIL ENERGIES IS FLAT
! *** FROM ZERO UP TO THE MAXIMUM 
      EREMAX=EREMAX*RN1
      ENEAFT=ENE-EREMAX
      ENDIF

      IF(ENE.LT.1.E-7)THEN
      CALL RANDM(RN1)
      ANGL=ACOS((2.*RN1)-1.)
      ELSE
      ETA=ACOS(SQRT((EREMAX/ENE)*(((1.+ANUMBE)**2.)/(4.*ANUMBE))))
      ANGL=ATAN((SIN(2.*ETA))/((1./ANUMBE)-COS(2.*ETA)))
      ENDIF
      
      RETURN
      END
! ************************************************
      SUBROUTINE ENERGY(ENE)
      COMMON IX
! *** USE VON NEUMANN'S METHOD ("HIT AND MISS") TO CALCULATE THE 
! *** ENERGY OF THE NEUTRON EMITED FROM A FISSION SOURCE.
   7  CALL RANDM(RN1)
! *** P=0 ABOVE ~10 MEV
      RN1=RN1*9.999
! *** APPLY MAXWELL SPECTRUM
      PROB=SQRT(RN1)*EXP(-RN1/1.4)
      CALL RANDM(RN2)
! *** THE MAXIMUM OF THE DISTRIBUTION AS WRITTEN IS LESS THAN 0.5 
      RN2=RN2*0.5
      IF(RN2.GT.PROB)GOTO 7    
      ENE=RN1   
      RETURN
      END
! ************************************************
      SUBROUTINE RANDM(RN1)
      COMMON IX
      IY=IX*6539
      IF(IY)5,6,6
    5 IY=IY+2147483647+1
    6 RN1=IY
      RN1=RN1*.4656613E-9
      IX=IY
      RETURN
      END
! ************************************************
      SUBROUTINE XSECT(ENE,XELAST,XABSOR)
! *** RETURNS THE APPROXIMATE CROSS SECTIONS (BARNS) FOR ABSORPTION AND 
! *** ELASTIC SCATTERING IN CARBON.
! *** HERE WE NEED THE ENERGY IN EV, NOT MEV
      X=ENE*1.E6

! *** CHECK THESE APPROXIMATIONS VS THE PLOTS IN THE HANDOUTS
      IF(X.LT.1.E4)THEN
      XELAST=5.
      ELSE
      XELAST=10.5-(1.346*LOG10(X))
      ENDIF
      XELAST=XELAST*1.E-24

      IF(X.LT.1.E3)THEN
      XABSOR=(6.442E-4)*(X**(-0.49421))
      ELSE IF(X.LT.1.E5) THEN
      XABSOR=1.5E-5
      ELSE IF(X.LT.5.E6) THEN
      XABSOR=2.E-5
      ELSE
      XABSOR=(4.E-06)*EXP(X*3.2189E-07) 
      ENDIF
      XABSOR=XABSOR*1.E-24
      RETURN
      END
! *************************************************
      SUBROUTINE EULER(EX,EY,EZ,ANGL,SX,SY,SZ)
!     THIS SUBROUTINE TAKES THE ORIGINAL LINEAR TRAJECTORY,
!     ROTATES IT TO LIE ALONG THE Z-AXIS, GENERATES A VECTOR
!     AT ZENITH ANGLE THETA = SCATTERING ANGLE 
!     AND AZIMUTHAL ANGLE FI = RANDOM * 2PI. THE
!     ORIGINAL AXIS IS NOW ROTATED BACK TAKING THE SCATTERING
!     VECTOR WITH IT. NOW WE HAVE THE SCATTERED
!     DIRECTION VECTOR (SX,SY,SZ).
!     WE USE EULER ANGLES TO PERFORM THE TRANSFORMATION.
      COMMON IX
! *** NORMALIZE THE DIRECTION TO A UNIT VECTOR (IN CASE IT WASN'T)
      S=SQRT(EX**2+EY**2+EZ**2)
      EX=EX/S
      EY=EY/S
      EZ=EZ/S
      BET=ACOS(EZ)
! *** THESE APPROXIMATIONS ARE ONLY NEEDED FOR COMPTON SCATTERING 
! *** FOR GAMMAS (BUT THEY WILL NOT HURT HERE)
      IF(ABS(BET).LT.0.027)ALF=0.0
      IF(ABS(BET).LT.0.027)GO TO 44
      ARG=EY/SIN(BET)
      AARG=ABS(ARG)
      IF(AARG.LT.1.0)GOTO 344
      ARG=ARG/(1.0001*AARG)
344   ALF=ASIN(ARG)
 44   CONTINUE
      SCO1=COS(ALF)*SIN(BET)+EX
      SCO1=ABS(SCO1)
      SCO2=ABS(EX)
      IF(SCO1.LT.SCO2)BET=-BET
      IF(SCO1.LT.SCO2)ALF=-ALF
      GAM=0.0
!     WE NOW HAVE THE EULER ANGLES OF ROTATION BETWEEN
!     Z-AXIS TO DIRECTION OF INITIAL PARTICLE.
      THET = ANGL
      CALL RANDM(RN1)
      FI = 6.2831853 * RN1
!     WE NOW HAVE SCATTERED THE PARTICLE FROM THE
!     Z-AXIS AND MUST ROTATE IT TO THE ORIGINAL UNSCATTERED
!     PARTICLE DIRECTION. CACULATE THE ROTATION MATRIX.
      R11 = COS(ALF)*COS(BET)*COS(GAM)-SIN(ALF)*SIN(GAM)
      R12 = COS(BET)*SIN(ALF)*COS(GAM)+COS(ALF)*SIN(GAM)
      R13 =-SIN(BET)*COS(GAM)
      R21 =-SIN(GAM)*COS(BET)*COS(ALF)-SIN(ALF)*COS(GAM)
      R22 =-SIN(GAM)*COS(BET)*SIN(ALF)+COS(ALF)*COS(GAM)
      R23 = SIN(BET)*SIN(GAM)
      R31 = SIN(BET)*COS(ALF)
      R32 = SIN(ALF)*SIN(BET)
      R33 = COS(BET)
      SOX = SIN(THET)*COS(FI)
      SOY = SIN(THET)*SIN(FI)
      SOZ = COS(THET)
      SX =  R11*SOX+R21*SOY+R31*SOZ
      SY =  R12*SOX+R22*SOY+R32*SOZ
      SZ =  R13*SOX+R23*SOY+R33*SOZ
!     WE NOW HAVE THE UNIT PROPAGATION VECTOR OF THE
!     SCATTERED PARTICLE IN THE *ORIGINAL* FRAME.
      RETURN
      END
! ************************************************

