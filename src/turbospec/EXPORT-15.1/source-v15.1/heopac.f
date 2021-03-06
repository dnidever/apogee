C
      SUBROUTINE HEOPAC(OMEGA,T,PROPAC)
C
      exp10(x)=exp(2.302585*x)
      A=7.02391+1.3380*ALOG10(T)
      A=exp10(A)
      A=1./A
      B=91.67+0.1033*T
      C=(15.57906-2.06158*ALOG10(T)-0.477352*(ALOG10(T))**2)/1.E7
      D=2.31317+3.8856E-4*T
      D=exp10(D)
      OMEGAC=274.3+.2762*T
      IF (OMEGA.GE.OMEGAC) GOTO 1
      OMEGAK=A*OMEGA**2*EXP(-OMEGA/B)
      GOTO 2
    1 OMEGAK=C*EXP(-OMEGA/D)
    2 CONTINUE
      OMEGAK=1.78*OMEGAK
      E=4.2432E-6-2.8854E-7*ALOG(T)
      F=1.2171E+5+258.28*T
      G=2.5830E-4-4.3429E-8*T
      H=1.1332E-2-1.1943E-3*ALOG(T)
      OMEGAP=-2973.3+600.73*ALOG(T)
      OMEGAT=1.5*OMEGAP
      IF(OMEGA.GT.OMEGAT) GOTO 13
      OMEGAX=E*EXP(-(OMEGA-OMEGAP)**2/F)
      GOTO 12
   13 OMEGAX=G*EXP(-H*OMEGA)
   12 CONTINUE
      OMEGAX=.10*OMEGAX
      DEL2=-4.033E+4+263.93*T
      DEL=SQRT(DEL2)
      APRIM=-7.7245+0.4246*ALOG10(T)
      APRIM=exp10(APRIM)
      AA=1./(1.125E+9+1.5866E+4*T+24.267*T**2)
      BB=1.2044+.4956*ALOG10(T)
      BB=exp10(BB)
      OMEGAZ=4161.1
      IF(OMEGA.GE.OMEGAZ) GOTO 22
      OMEGAY=(APRIM*DEL*OMEGA*EXP((OMEGA-OMEGAZ)/(.6952*T)))/((OMEGA-
     & OMEGAZ)**2 +DEL2)
      GOTO 23
   22 IF(OMEGA.GT.(OMEGAZ+1.5*DEL)) GOTO 24
      OMEGAY=(APRIM*DEL*OMEGA)/((OMEGA-OMEGAZ)**2+DEL2)
      GOTO 23
   24 OMEGAY=AA*OMEGA*EXP(-(OMEGA-OMEGAZ)/BB)
   23 CONTINUE
      PROPAC=OMEGAK+OMEGAX+OMEGAY
      RETURN
      END
