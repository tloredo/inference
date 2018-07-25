C This version has the kmdiff call within kmestm.

C      **********************************************************************
C      ********************* SUBROUTINE KMESTM ******************************
C      **********************************************************************
C
       SUBROUTINE KMESTM(IND,X,NTOT, IERR,SX,VX,SMEAN,ERROR,
     +                   KDIFF,START,BINSIZ,LSTEP,
     +                   IXU,ZXU,IXC,ZXC,BWRK1,BWRK2,BWRK3,
     +                   Z1,ITEMP,INDEXX,ATIE,RISK,IWRK1,WRK1)
C
C      *       THIS SUBROUTINE COMPUTES THE PL ESTIMATOR, MEAN AND ITS      *
C      *       ERROR FOR THE X VARIABLE.                                    *
C      *                                                                    *
C      *       INPUT IND(I): INDICATOR OF CENSORING                         *
C      *               X(I): DATA                                           *
C      *                NTOT : NO. OF DATA POINTS                           *
C      *               KDIFF : IF KDIFF = 1, PRINT OUT DIFFERENTIAL KM      *
C      *               START : STARTING POINT OF BINING                     *
C      *               BINSIZ: WIDTH OF THE BIN                             *
C      *               LSTEP : NUMBER OF BINS                               *
C      *              ATIE(I): NUMBER OF TIED DATA POINTS AT ITH DATA VALUE *
C      *              RISK(I): RISK SET FOR THE ITH DATA VALUE              *
C      *                                                                    *
C      *       WORK      ZXU : DATA ARRAY CONTAINING THE UNCENSORED POINTS  *
C      *                 ZXC : DATA ARRAY CONTAINING THE CENSORED POINTS    *
C      *                 IXU : NO. OF UNCENSORED DATA POINTS                *
C      *                 IXC : NO. OF CENSORED DATA POINTS                  *
C      *              ICHANGE: IF THE LAST VALUE IS CENSORED, THE VALUE     *
C      *                       NEEDS TO BE CHANGED TO A DETECTION.          *
C      *                       THIS INDICATOR IS SET TO BE -1,IF THE LAST   *
C      *                       VALUE IS CHANGED.                            *
C      *                                                                    *
C      *       OUTPUT    SX  : PL ESTIMATOR                                 *
C      *                 VX  : ERROR OF PL ESTIMATOR                        *
C      *                SMEAN: MEAN                                         *
C      *                ERROR: ERROR OF THE MEAN                            *
C      *                                                                    *
C      *     SUBROUTINES                                                    *
C      *                 SORT1, XVAR, PLESTM, KMADJ, KMPRNT                 *
C      *                                                                    *
C
C       NOTE:  ZXU and ZXC are dim(NTOT), but contain only IXU and IXC items
C       on output.  The Python caller should adjust them.  SUMRY (below)
C       expects them to be of size IXU and IXC.

c Basic KM params & results:
cf2py intent(in) NTOT, IND, X
cf2py intent(out) IERR, SX, VX, SMEAN, ERROR
cf2py intent(out) IXU, ZXU, IXC, ZXC

c Diff'l KM params:
cf2py intent(in) KDIFF, START, BINSIZ, LSTEP
cf2py intent(out) BWRK1, BWRK2, BWRK3

c Workspace:
cf2py intent(cache,hide) Z1,ITEMP,INDEXX,ATIE,RISK,IWRK1,WRK1

C       For diff'l KM, BWRK1-3 are the le
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION IND(NTOT),X(NTOT),ZXU(NTOT),ZXC(NTOT)
       DIMENSION SX(NTOT),VX(NTOT),Z1(NTOT),ITEMP(NTOT)
       DIMENSION INDEXX(NTOT),ATIE(NTOT),RISK(NTOT)
       DIMENSION BWRK1(LSTEP),BWRK2(LSTEP),BWRK3(LSTEP)
       DIMENSION IWRK1(NTOT),WRK1(NTOT)

C
C      *       DISTINGUISH UNCENSORED AND CENSORED DATA AND SORT THEM IN    *
C      *       ASCENDING ORDER. THEN CALL PL ESTIMATOR SUBROUTINE           *
C      *       "PL".                                                        *
C

C
C      *    XVAR DISTINGUISHES DETECTED POINTS AND CENSORED POINTS           *
C
       CALL XVAR(IND,X,NTOT,ISIGN,ZXU,ZXC,IXU,IXC,IWRK1,
     +           ATIE,RISK,WRK1,Z1,LTOT,INDEXX,IERR)
C
       IF (IERR .EQ. -1) RETURN

cc          IF(ISIGN.EQ.1) PRINT 56,NTOT,IXC
cc          IF(ISIGN.EQ.-1) PRINT 57,NTOT,IXC
cc          PRINT *
cc   56  FORMAT(8X,'# OF DATA POINTS :',I4,' # OF LOWER LIMITS :',I4)
cc   57  FORMAT(8X,'# OF DATA POINTS :',I4,' # OF UPPER LIMITS :',I4)

C
C      *  SET A FEW DUMMY ARRAYS TO USE THE SORTING PRGRAM                  *
C
       ISTEMP = ISIGN
   59  DO 60 I=1,NTOT
          ITEMP(I)=0
          Z1(I)=1.0
   60  CONTINUE
C
       CALL SORT1(ITEMP,Z1,ZXU,IXU,INDEXX)
C
       CALL SORT1(ITEMP,Z1,ZXC,IXC,INDEXX)
C
C      *  CALL SUBROUTINE "PLESTM" TO COMPUTE KM ESTIMATOR                  *
C
       CALL PLEST1(ZXU,ZXC,IXU,IXC,SX,VX,NTOT,SMEAN,ERROR,ICHANGE,
     +             NCHANGE,IWRK1)
C
C      *        IF THE DATA CONTAINS UPPER LIMITS, CHANGE THE               *
C      *        SIGN OF THE MEAN.                                           *
C
       ISIGN = ISTEMP
       SMEAN=ISIGN*SMEAN
C
C      * SUBROUTINE KMADJ IS CALLED TO ADJUST THE PRODUCT-LIMIT ESTIMATOR   *
C      * BACK TO THE ORIGIONAL CENSORING PATTERN AND TO REMOVE TIES.        *
C
       CALL KMADJ(ZXU,ZXC,NTOT,IXU,IXC,SX,ISIGN,NTEMP,WRK1,VX)

cc       CALL KMPRNT(ZXU,ZXC,NTOT,NTEMP,IXU,IXC,SX,VX,ISIGN,OUTPUT,
cc     +             ICHANGE,SMEAN,ERROR,IPRINT)


C      *  SUBROUTINE KMDIF IS CALLED IF THE USER HAS REQUESTED A            *
C      *  DIFFERENTIAL KM ESTIMATOR.                                        *

       IF(KDIFF .EQ. 1) THEN

          CALL KMDIF(SX,ZXU,BWRK1,BWRK2,BWRK3,WRK1,NTOT,START,
     +               BINSIZ,LSTEP,IXU)

       ENDIF

C
C      *   IF THE LAST VALUE WAS CHANGED FROM AN UPPER LIMIT TO A          *
C      *   DETECTION, CHANGE THE NUMBER BACK TO ITS ORIGINAL VALUE.        *
C
       IF(ICHANGE.EQ.-1) THEN
          IXU=IXU-NCHANGE
          IXC=IXC+NCHANGE
       ENDIF
       
       RETURN
       END

C
C      **********************************************************************
C      ******************** SUBROUTINE KMADJ  *******************************
C      **********************************************************************
C
       SUBROUTINE KMADJ(ZU,ZC,NTOT,IU,IC,S,ISIGN,NTEMP,F,V)
C
C      *       THIS SUBROUTINE RESTORES THE DATA AND THE PRODUCT-LIMIT      *
C      *       ESTIMATOR TO UPPER-LIMITS FORM, IF THE DATA WAS IN THAT FORM *
C      *       INITIALLY.  TIES AT CENSORED POINTS ARE ALSO REMOVED TO      *
C      *       MAKE THE PRINTOUT CONSISTENT.                                *
C      *                                                                    *
C      *       INPUT    ZU(I)  :  UNCENSORED DATA POINTS                    *
C      *                ZC(I)  :  CENSORED DATA POINTS                      *
C      *                NTOT   :  TOTAL NUMBER OF DATA POINTS,=IU+IC.       *
C      *                 IU    :  NUMBER OF UNCENSORED DATA POINTS          *
C      *                 IC    :  NUMBER OF CENSORED DATA POINTS            *
C      *                 S(L)  :  PL ESTIMATOR                              *
C      *       OUTPUT  NTEMP   :  VALUE OF NTOT AFTER ADJUSTMENT FOR TIES   *
C      *       OTHER    F      :  PROBABILITY MASS ASSIGNED TO EACH         *
C      *                             UNCENSORED POINT(=JUMP IN S AT THE     *
C      *                                                  POINT)            *
C      *                                                                    *

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION ZU(NTOT),ZC(NTOT),S(NTOT),F(NTOT),V(NTOT)

C
C      *  FOR LEFT-CENSORED DATASETS (I.E. UPPER LIMITS),                  *
C      *  CALCULATE JUMP IN SURVIVOR FUNCTION AT UNCENSORED POINTS         *
C
       IF(ISIGN.LT.0) THEN
            F(1) = 1.0 - S(1)
            DO 120 I = 2, IU
               F(I) = S(I-1)-S(I)
120         CONTINUE
C
C      *  REVERSE SEQUENCE OF UNCENSORED POINTS, JUMPS AND ERRORS         *
C
            J = IU/2
            DO 150 I =1, J

               Z = ZU(I)*(-1.0)
               ZU(I) = ZU(IU-I+1)*(-1.0)
               ZU(IU-I+1) = Z

               FTEMP = F(I)
               F(I) = F(IU-I+1)
               F(IU-I+1) = FTEMP 

               VTEMP = V(I)
C               V(I) = V(IU-I+1)
C               V(IU-I+1) = VTEMP
               V(I) = V(IU-I)
               V(IU-I) = VTEMP

150         CONTINUE

            IF(2*J.NE.IU) THEN
               ZU(J+1) = ZU(J+1)*(-1.0)
            ENDIF

C
C      *  REVERSE SEQUENCE OF CENSORED POINTS                              *
C
            J = IC/2
            DO 155 I = 1, J
               Z = ZC(I) * (-1.0)
               ZC(I) = ZC(IC-I+1)*(-1.0)
               ZC(IC-I+1) = Z
155         CONTINUE
 
            IF(2*J.NE.IC) THEN
               ZC(J+1) = ZC(J+1)*(-1.0)
            ENDIF

C
C      *  COMPUTE SURVIVOR FUNCTION FOR REVERSED DATA                     *
C
            DO 170 I = 1, IU
               S(I) = 1
               DO 160 J = 1, I
                  S(I) = S(I) - F(J)
160            CONTINUE
170         CONTINUE
         ENDIF   

C      *   CORRECTS FOR TIES AT THE UNCENSORED POINTS                      *
C      *   NOTICE THAT IF THERE ARE TIES AT THESE POINTS, THEN BOTH        *
C      *   IU AND NTEMP ARE REDUCED.                                       *

       NTEMP = NTOT
       K = 1
190    IF(ZU(K).EQ.ZU(K+1)) THEN
          DO 195 I = K, IU-1
             ZU(I)=ZU(I+1)
             S(I)=S(I+1)
             V(I) = V(I+1)               
195       CONTINUE
          IU = IU -1
          NTEMP = NTEMP - 1
       ELSE
          K  = K + 1
       ENDIF
       IF(K.LT.IU) GOTO 190

       RETURN
       END


C
C      **********************************************************************
C      ******************** SUBROUTINE KMDIF ********************************
C      **********************************************************************
C
       SUBROUTINE KMDIF(S,ZU,BS,BL,DIFF,F,NTOT,START,BINSIZ,LSTEP,IU)

C
C      *       THIS SUBROUTINE COMPUTES AND PRINTS THE DIFFERENTIAL KM      *
C      *       ESTIMATOR BASED ON WARDLE AND KNAPP, 'THE STATISTICAL        *
C      *       DISTRIBUTION OF THE NEUTRAL-HYDROGEN CONTENT OF S0 GALAXIES',*
C      *       ASTRN. J., 91:23 1986.                                       *
C      *                                                                    *
C      *       INPUT    ZU(I)  :  UNCENSORED DATA POINTS                    *
C      *                NTOT   :  TOTAL NUMBER OF DATA POINTS,=IU+IC.       *
C      *                 IU    :  NUMBER OF UNCENSORED DATA POINTS          *
C      *                 S(L)  :  PL ESTIMATOR                              *
C      *                 OUT   :  OUTPUT FILE NAME. IF IT IS BLANK          *
C      *                          THE RESULTS WILL BE SHOWN ON THE          *
C      *                          TERMINAL.                                 *
C      *                START  :  STARTING VALUE OF THE FIRST BIN           *
C      *                BINSIZ :  WIDTH OF THE BIN                          *
C      *                LSTEP  :  NUMBER OF BINS                            *
C      *              ICHANGE  :  INDICATES IF THE LAST POINT (OR THE       *
C      *                            FIRST POINT FOR UPPER LIMITS DATA)      *
C      *                            HAS BEEN CHANGED TO A DETECTION.        *
C      *                                                                    *
C      *      OTHERS                                                        *
C      *               BS(J)   :  STARTING VALUE FOR THE BIN J              *
C      *               BL(J)   :  ENDING VALUE FOR THE BIN J                *
C      *               DIFF(J) :  DIFFERENTIAL KM ESTIMATOR AT BIN J        *
C      *               F(I)    :  MASS OF THE I TH DATA POINT               *
C      *                                                                    *
C      *                                                                    *
C
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION ZU(NTOT),S(NTOT),F(NTOT)
       DIMENSION BS(LSTEP),BL(LSTEP),DIFF(LSTEP)
C
C       *    FIRST, COMPUTE THE MASS OF EACH POINT                          *
C
        F(1) = 1.0 -S(1)
       
        DO 610 I = 2, IU
           F(I) = DABS(S(I) - S(I-1))
  610   CONTINUE

C
C       *  SET THE BIN BOUNDARIES.                                          *
C
        DO 620 J = 1, LSTEP
           BS(J) = START + BINSIZ*(J-1)
           BL(J) = START + BINSIZ*J
  620   CONTINUE

        I = 1
        J = 0

  630   J = J + 1
        DIFF(J) = 0.0

C
C      *       CHECK WHETHER THE I-TH POINT IS IN THE BIN                  *
C
  640  IF(J .LE. LSTEP) THEN
         IF(ZU(I) .LT. BS(J)) THEN
            IF(I .GE. IU) THEN
               GOTO 630
            ENDIF
            I = I + 1
            GOTO 640
         ENDIF

C      *       COMPUTE THE DIFFERENTIAL KM                                 *
C
         IF(ZU(I) .GE. BL(J)) GOTO 630
         DIFF(J) = DIFF(J) + F(I)
   
         IF(I .LT. IU) THEN
            I = I + 1
            GOTO 640
         ENDIF
         GOTO 630
       ENDIF

C
C      * MULTIPLY DIFF(I) BY THE TOTAL NUMBER OF POINTS TO GET THE NUMBER    *
C      * OF POINTS IN EACH BIN.                                              *
C
          SUM = 0.0
          DO 690 I = 1, LSTEP
             DIFF(I) =DIFF(I)*NTOT
             CENTER = 0.5*(BS(I) + BL(I))
             SUM = SUM + DIFF(I)
  690     CONTINUE

cc             PRINT 700, SUM

       RETURN
       END


C
C      **********************************************************************
C      ********************* SUBROUTINE PLESTM  *****************************
C      **********************************************************************
C
       SUBROUTINE PLESTM(NU,U,NC,C,S,V,SMEAN,SIGMA,ICHANGE,
     +                   NCHANGE,L)
C       
C      *      THIS SUBROUTINE COMPUTES PL ESTIMATOR AND THE MEAN            *
C      *      AND ITS ERROR.                                                *
C      *                                                                    *
C      *       INPUT     U : UNCENSORED DATA POINTS                         *
C      *                 C : CENSORED DATA POINTS                           *
C      *                NU : NO. OF UNCENSORED DATA POINTS                  *
C      *                NC : NO. OF CENSORED DATA POINTS                    *
C      *               NTOT: TOTAL NUMBER OF DATA POINTS                    *
C      *                                                                    *
C      *       WORK      L : RANK OF THE UNCENSORED DATA                    *
C      *               VAR : VARIANCE OF THE MEAN                           *
C      *                KD : NUMBER OF TIED DATA POINTS                     *
C      *                                                                    *
C      *       OUTPUT    S : PL ESTIMATOR                                   *
C      *                 V : ERROR FOR THE PL ESTIMATOR                     *
C      *             SMEAN : MEAN OF THE DATA                               *
C      *             SIGMA : ERROR OF THE MEAN                              *
C      *            ICHANGE: IF THE LAST VALUE IS CENSORED, WE NEED TO      *
C      *                     CHANGE IT TO A DETECTION. THEN ICHANGE=-1,     *
C      *                     OTHERWISE ICHANGE=1.                           *
C      *            NCHANGE: IF ICHANGE = -1 AND THE LAST VALUE IS TIED     *
C      *                     WITH OTHER CENSORED VALUES, THIS RECORDS THE   *
C      *                     NUMBER OF TIED OBSERVATIONS (ALL OF THEM NEED  *
C      *                     TO BE CHANGED TO DETECTIONS).                  *
C      *                                                                    *
C      *       FIRST HALF OF THE PROGRAM IS FROM ELISA T. LEE, "STATISTICAL *
C      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING *
C      *       PUBLICATIONS (BELMONT:CA); WITH THE GRAPHIC ROUTINES REMOVED.*
C      *       FORMULAS USED FOR COMPUTATION OF THE MEAN AND THE ERROR ARE  *
C      *       FROM RUPERT G. MILLER, "SURVIVAL ANALYSIS", 1981,            *
C      *       JOHN WILEY & SONS (NY:NY)                                    *
C      *                                                                    *
C

cf2py intent(in) NU, U, NC, C
cf2py intent(out) S, V, SMEAN, SIGMA, ICHANGE, NCHANGE
cf2py intent(cache,hide) L

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION U(NU),C(NC),S(NU),V(NU),L(NU)
C
C      *          COMPUTE THE RANK (L) OF THE UNCENSORED POINTS             *
C
C*******     IF THE LAST VALUE IS CENSORED, CHANGE IT TO A DETECTION        *
C

C      THE FOLLOWING LOOP HAS BEEN MODIFIED AND NCHANGE ADDED TO THE 
C      PROGRAM TO COVER THE CASE WHEN TIED NONDETECTIONS ARE THE LARGEST
C      VALUE.  MODIFIED 4/92

       ICHANGE=1
       NCHANGE = 0
 13    IF(NC .NE. 0)THEN 
          IF(U(NU) .LE. C(NC))THEN 
             U(NU+1)=C(NC)
             NU=NU+1
             NC=NC-1
             NCHANGE = NCHANGE + 1
             ICHANGE=-1
          ELSE
             GOTO 15
          ENDIF
       ELSE
          GOTO 15
       ENDIF
       GOTO 13
C
 15    K=1
       KK=0
       NT=NU+NC
       IF(NC .NE. 0) THEN 
          DO 10 I=1,NU
             IF(KK .NE. NC) THEN
                DO 20 J=K,NC
                   K1=J
                   IF(C(J) .GE. U(I)) GOTO 1
                   KK=KK+1
  20            CONTINUE
             ENDIF
   1         K=K1
             L(I)=I+KK
  10      CONTINUE
       ELSE
          DO 19 I=1,NU
             L(I)=I
  19      CONTINUE
       ENDIF
C
C      *       COMPUTE P(T) FOR ALL UNCENSORED POINTS BASED ON RANK (L)     *
C
       V1=0.0
       PT=1.0
       XNT=NT
       DO 12 I=1,NU
          XL=L(I)
          PT=PT*((XNT-XL)/(XNT-XL+1.0))
          S(I)=PT
          IF((XNT-XL) .LE. 0.0) THEN
             VV=0.0
          ELSE
             V1=V1+1.0/((XNT-XL)*(XNT-XL+1.0))
             VV=DSQRT((PT**2)*V1)
          ENDIF
          V(I)=VV
  12   CONTINUE

C
C      *        COMPUTE THE MEAN                                            *
C      *        REF. FOR THE MEAN AND ERROR :                               *
C      *          MILLER, R. G. JR. 1981, "SURVIVAL ANALYSIS"               *
C      *          PP. 70-71 AND 195-198.                                    *
C
       SMEAN=U(1)
       I=2
  30   K=0
  40   IF((U(I+K).NE.U(I-1)).AND.(I+K.LE.NU)) THEN
          SMEAN=SMEAN+S(I+K-1)*(U(I+K)-U(I-1))
          IF(I+K.LT.NU) THEN
             I=I+K+1
             GOTO 30
          ENDIF
       ELSEIF(U(I+K).EQ.U(I-1)) THEN
          K=K+1
          GOTO 40
       ENDIF
C
C      *              COMPUTE THE ERROR OF THE MEAN                         *
C
       J=2    
       VAR=0.0
  120  I=J
       SSUM=0
  130  K=0
  140  IF((U(I+K).EQ.U(I-1)).AND.(I+K.LE.NU)) GOTO 145
          IF(U(I+K).EQ.U(I-1)) THEN
             K=K+1
             GOTO 140
          ENDIF
  145     SSUM=SSUM+S(I+K-1)*(U(I+K)-U(I-1))
          IF(I+K.LT.NU) THEN
             I=I+K+1
             GOTO 130
          ENDIF
C
C      *          KD IS NO. OF TIED OBSERVATIONS AT THAT POINT              *
C
       KD=1
  180  IF(U(J-1+KD).LE.U(J-1)) THEN
          KD=KD+1
          GOTO 180
       ENDIF
       XL=L(J-1)
       D=KD
       B=XNT-XL-D+1.0
C
C      *       IF THE LAST FEW DATA POINTS ARE UNCENSORED AND TIED, SKIP    *
C      *       THE NEXT LINES TO AVOID DIVISION BY 0.                       *
C
       IF(B .NE. 0.0) THEN
          VAR=VAR+SSUM*SSUM*D/((XNT-XL+1)*B)
          J=J+KD
          IF(J.LE.NU) GOTO 120
       ENDIF
  200  SIGMA=DSQRT(VAR)
   
       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE PLEST1  *****************************
C      **********************************************************************
C
C      This is PLESTM from the original ASURV; "1" refers to its use of a
C      single array dimension for the input args.  PLESTM above is the same
C      subroutine, but with the input args having individual dimensions, to
C      make direct calling from Python more natural.
C      **********************************************************************
C
       SUBROUTINE PLEST1(U,C,NU,NC,S,V,NTOT,SMEAN,SIGMA,ICHANGE,
     +                   NCHANGE,L)
C       
C      *      THIS SUBROUTINE COMPUTES PL ESTIMATOR AND THE MEAN            *
C      *      AND ITS ERROR.                                                *
C      *                                                                    *
C      *       INPUT     U : UNCENSORED DATA POINTS                         *
C      *                 C : CENSORED DATA POINTS                           *
C      *                NU : NO. OF UNCENSORED DATA POINTS                  *
C      *                NC : NO. OF CENSORED DATA POINTS                    *
C      *               NTOT: TOTAL NUMBER OF DATA POINTS                    *
C      *                                                                    *
C      *       WORK      L : RANK OF THE UNCENSORED DATA                    *
C      *               VAR : VARIANCE OF THE MEAN                           *
C      *                KD : NUMBER OF TIED DATA POINTS                     *
C      *                                                                    *
C      *       OUTPUT    S : PL ESTIMATOR                                   *
C      *                 V : ERROR FOR THE PL ESTIMATOR                     *
C      *             SMEAN : MEAN OF THE DATA                               *
C      *             SIGMA : ERROR OF THE MEAN                              *
C      *            ICHANGE: IF THE LAST VALUE IS CENSORED, WE NEED TO      *
C      *                     CHANGE IT TO A DETECTION. THEN ICHANGE=-1,     *
C      *                     OTHERWISE ICHANGE=1.                           *
C      *            NCHANGE: IF ICHANGE = -1 AND THE LAST VALUE IS TIED     *
C      *                     WITH OTHER CENSORED VALUES, THIS RECORDS THE   *
C      *                     NUMBER OF TIED OBSERVATIONS (ALL OF THEM NEED  *
C      *                     TO BE CHANGED TO DETECTIONS).                  *
C      *                                                                    *
C      *       FIRST HALF OF THE PROGRAM IS FROM ELISA T. LEE, "STATISTICAL *
C      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING *
C      *       PUBLICATIONS (BELMONT:CA); WITH THE GRAPHIC ROUTINES REMOVED.*
C      *       FORMULAS USED FOR COMPUTATION OF THE MEAN AND THE ERROR ARE  *
C      *       FROM RUPERT G. MILLER, "SURVIVAL ANALYSIS", 1981,            *
C      *       JOHN WILEY & SONS (NY:NY)                                    *
C      *                                                                    *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION U(NTOT),C(NTOT),S(NTOT),V(NTOT),L(NTOT)
C
C      *          COMPUTE THE RANK (L) OF THE UNCENSORED POINTS             *
C
C*******     IF THE LAST VALUE IS CENSORED, CHANGE IT TO A DETECTION        *
C

C      THE FOLLOWING LOOP HAS BEEN MODIFIED AND NCHANGE ADDED TO THE 
C      PROGRAM TO COVER THE CASE WHEN TIED NONDETECTIONS ARE THE LARGEST
C      VALUE.  MODIFIED 4/92

       ICHANGE=1
       NCHANGE = 0
 13    IF(NC .NE. 0)THEN 
          IF(U(NU) .LE. C(NC))THEN 
             U(NU+1)=C(NC)
             NU=NU+1
             NC=NC-1
             NCHANGE = NCHANGE + 1
             ICHANGE=-1
          ELSE
             GOTO 15
          ENDIF
       ELSE
          GOTO 15
       ENDIF
       GOTO 13
C
 15    K=1
       KK=0
       NT=NU+NC
       IF(NC .NE. 0) THEN 
          DO 10 I=1,NU
             IF(KK .NE. NC) THEN
                DO 20 J=K,NC
                   K1=J
                   IF(C(J) .GE. U(I)) GOTO 1
                   KK=KK+1
  20            CONTINUE
             ENDIF
   1         K=K1
             L(I)=I+KK
  10      CONTINUE
       ELSE
          DO 19 I=1,NU
             L(I)=I
  19      CONTINUE
       ENDIF
C
C      *       COMPUTE P(T) FOR ALL UNCENSORED POINTS BASED ON RANK (L)     *
C
       V1=0.0
       PT=1.0
       XNT=NT
       DO 12 I=1,NU
          XL=L(I)
          PT=PT*((XNT-XL)/(XNT-XL+1.0))
          S(I)=PT
          IF((XNT-XL) .LE. 0.0) THEN
             VV=0.0
          ELSE
             V1=V1+1.0/((XNT-XL)*(XNT-XL+1.0))
             VV=DSQRT((PT**2)*V1)
          ENDIF
          V(I)=VV
  12   CONTINUE

C
C      *        COMPUTE THE MEAN                                            *
C      *        REF. FOR THE MEAN AND ERROR :                               *
C      *          MILLER, R. G. JR. 1981, "SURVIVAL ANALYSIS"               *
C      *          PP. 70-71 AND 195-198.                                    *
C
       SMEAN=U(1)
       I=2
  30   K=0
  40   IF((U(I+K).NE.U(I-1)).AND.(I+K.LE.NU)) THEN
          SMEAN=SMEAN+S(I+K-1)*(U(I+K)-U(I-1))
          IF(I+K.LT.NU) THEN
             I=I+K+1
             GOTO 30
          ENDIF
       ELSEIF(U(I+K).EQ.U(I-1)) THEN
          K=K+1
          GOTO 40
       ENDIF
C
C      *              COMPUTE THE ERROR OF THE MEAN                         *
C
       J=2    
       VAR=0.0
  120  I=J
       SSUM=0
  130  K=0
  140  IF((U(I+K).EQ.U(I-1)).AND.(I+K.LE.NU)) GOTO 145
          IF(U(I+K).EQ.U(I-1)) THEN
             K=K+1
             GOTO 140
          ENDIF
  145     SSUM=SSUM+S(I+K-1)*(U(I+K)-U(I-1))
          IF(I+K.LT.NU) THEN
             I=I+K+1
             GOTO 130
          ENDIF
C
C      *          KD IS NO. OF TIED OBSERVATIONS AT THAT POINT              *
C
       KD=1
  180  IF(U(J-1+KD).LE.U(J-1)) THEN
          KD=KD+1
          GOTO 180
       ENDIF
       XL=L(J-1)
       D=KD
       B=XNT-XL-D+1.0
C
C      *       IF THE LAST FEW DATA POINTS ARE UNCENSORED AND TIED, SKIP    *
C      *       THE NEXT LINES TO AVOID DIVISION BY 0.                       *
C
       IF(B .NE. 0.0) THEN
          VAR=VAR+SSUM*SSUM*D/((XNT-XL+1)*B)
          J=J+KD
          IF(J.LE.NU) GOTO 120
       ENDIF
  200  SIGMA=DSQRT(VAR)
   
       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE QUART  ******************************
C      Altered from ASURV's SUMRY so U, S have dim(IU).
C      **********************************************************************
C
       SUBROUTINE QUART(U,S,IU,FINT)
C
C      *       THIS SUBROUTINE CALCULATES AND PRINTS 75, 50, AND            *
C      *       25, PERCENTILES  OF SURVIVAL FOR A SURVIVAL CURVE.           *
C      *       S AND U ARE ARRAYS CONTAINING POINTS FOR WHICH VALUES OF THE *
C      *       SURVIVAL CURVE WERE CALCULATED, IU IS THE NUMBER OF          *
C      *       UNCENSORED POINTS.  ADOPTED FROM ELISA T. LEE, "STATISTICAL  *
C      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME LEARNING *
C      *       PUBLICATIONS (BELMONT:CA).                                   *
C      *                                                                    *
C      *       INPUT       U : UNCENSORED DATA                              *
C      *                   S : PL ESTIMATOR                                 *
C      *       WORK       TY : PERCENTILE INDICATOR AT 75, 50, 25           *
C      *       OUTPUT    FINT: VALUES AT 75, 50, 25 PERCENTILES             *
C      *                                                                    *
C      *     THIS SUBROUTINE IS FROM SMSDA, EXCEPT PRINTING COMMAND.        *
C      *                                                                    *
C
C

cf2py intent(in) U, S, IU
cf2py intent(out) FINT

       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION U(IU),S(IU),TY(3),FINT(3)
C
C
       TY(1)=0.75
       TY(2)=0.50
       TY(3)=0.25

C      *       IF THE NUMBER OF DATA POINTS IS SMALLER THAN 4, NO          *
C      *       PERCENTILES CAN BE OBTAINED.                                *
C
       DO 40 I=1,3
          FINT(I)=0.0
  40   CONTINUE

       L=1
       IF(IU.LE.3) RETURN
       NN=IU-1

       DO 100 I=1,3
          IF(S(1).LT.TY(I)) THEN
             FINT(I) = U(1)-(TY(I)-S(1))/(1-S(1))*(U(1)-0)
          ELSE
             DO  90 J=L,NN
                IF((S(J).GE.TY(I)) .AND. (S(J+1).LE.TY(I))) THEN 
                   FINT(I)=U(J)-(S(J)-TY(I))/(S(J)-S(J+1))*(U(J)-U(J+1))
                   L=J+1
                   GOTO 100
                ENDIF
  90         CONTINUE
          ENDIF
 100   CONTINUE

       RETURN
       END

C
C      **********************************************************************
C      ********************* SUBROUTINE XVAR  *******************************
C      **********************************************************************
C
       SUBROUTINE XVAR(IND,X,NTOT,ISIGN,ZU,ZC,IU,IC,ISAVE,
     +                 ATIE,RISK,XT,ZTEMP,LTOT,INDEXX,IERR)
C
C      *       THIS SUBROUTINE DISTINGUISHES UNCENSORED AND CENSORED        *
C      *       DATA IN THE X VARIABLE AND SORTS IT INTO ZU AND ZC. ALSO,    *
C      *       IF THE DATA CONTAIN UPPER LIMITS, THE SIGN OF THE            *
C      *       VALUES ARE CHANGED SO THAT THE PROGRAM FOR THE LOWER         *
C      *       LIMITS CAN BE USED. ADOPTED FROM ELISA T. LEE, "STATISTICAL  *
C      *       METHODS FOR SURVIVAL DATA ANALYSIS", 1980, LIFETIME          *
C      *       LEARNING PUBLICATIONS (BELMONT:CA).                          *
C      *                                                                    *
C      *       INPUT      IND(I) : INDICATOR OF CENSORING                 *
C      *                    X(I) : VARIABLE                               *
C      *                   NTOT    : TOTAL NUMBER OF DATA POINTS            *
C      *                                                                    *
C      *       OUTPUT      ISIGN   : IF LOWER LIMIT, ISIGN = 1              *
C      *                             IF UPPER LIMIT, ISIGN = -1             *
C      *                   ZU(K)   : UNCENSORED DATA POINTS IN X(J,I)       *
C      *                   ZC(K)   : CENSORED DATA POINTS IN X(J,I)         *
C      *                    IU     : NUMBER OF UNCENSORED DATA POINTS       *
C      *                    IC     : NUMBER OF CENSORED DATA POINTS         *
C      *                   RISK(L) : RISK SET                               *
C      *                  ATIE(L)  : NUMBER OF TIED DATA POINTS             *
C      *                                                                    *
C      *       OTHER                                                        *
C      *                  ISAVE(I) : TEMPORARY SAVE OF ISIGN FOR EACH POINT *
C      *                             ALSO USED FOR TEMPORARY CENSORSHIP     *
C      *                             INDICATOR                              *
C      *                   XT(I)   : = X(J,I)                               *
C      *                                                                    *
C      *        SUBROUTINES                                                 *
C      *                   SORT1                                            *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

       DIMENSION IND(NTOT),X(NTOT),ZU(NTOT)
       DIMENSION ZC(NTOT),ISAVE(NTOT),ATIE(NTOT),RISK(NTOT)
       DIMENSION XT(NTOT),ZTEMP(NTOT)
       DIMENSION INDEXX(NTOT)
C
       ISIGN=1     
       IC=0
       IU=0
       IERR = 0
C
C      *    FIND THE CENSORSHIP OF THE DATA SET. -1 FOR THE UPPER LIMIT    *
C      *    AND 1 FOR THE LOWER LIMIT                                      *
C
       DO 100 I=1,NTOT
          ISAVE(I) = 0
          IF(IND(I) .EQ. 0) GOTO 100
          ISIGN=IND(I)/IABS(IND(I))
          ISAVE(I) = ISIGN
  100  CONTINUE

C      * CHECK WHETHER THE UPPER AND LOWER LIMITS ARE MIXED IN ONE        *
C      * VARIABLE. IF SO, THE PROGRAM IS TERMINATED.                      *
C
       DO 110 I = 1, NTOT
          IF(ISAVE(I) .EQ. 0) GOTO 110
          IF(ISAVE(I) .NE. ISIGN) THEN
             IERR = -1
             RETURN
cc             PRINT *,'YOU CANNOT HAVE BOTH UPPER AND LOWER LIMITS'
cc             PRINT *,'IN ONE VARIABLE AT THE SAME TIME.'
cc             PRINT *,'PLEASE CHECK THE DATA. THE PROGRAM IS TERMINATED.'
          ENDIF
  110  CONTINUE
C
C      *    IN CASE THE DATA HAS UPPER LIMITS IT IS MULTIPLIED BY ISIGN   *
C      *    TO MAKE THE DATA HAVE LOWER LIMITS (RIGHT CENSORSHIP).        *
C
       DO 280 L = 1, NTOT
          ATIE(L) = 0.0
          XT(L) = REAL(ISIGN)*X(L)
          ZTEMP(L) = 0.0
          ISAVE(L) = IND(L)
  280  CONTINUE
C
C     *     DATA POINTS ARE ARRANGED FROM SMALLEST TO LARGEST.             *
C     *     DETECTED AND CENSORED DATA POINTS ARE SEPARATED.               *
C     *     THEN RISK SETS AND TIED DATA POINTS ARE FOUND.                 *
C

       CALL SORT1(ISAVE,ZTEMP,XT,NTOT,INDEXX)

       L = 1

       DO 300 I=1,NTOT
          K=IABS(ISAVE(I))
          IF(K .EQ. 0) THEN 
              IU=IU+1
              ZU(IU)= XT(I)
              IF(IU .NE. 1) THEN
                 IF(ZU(IU) .EQ. ZU(IU-1)) THEN
                    ATIE(L) = ATIE(L) + 1.0
                    RISK(L) = REAL(NTOT - I)
                 ELSE
                    ATIE(L) = 1.0
                    RISK(L) = REAL(NTOT - I)
                    L = L + 1
                 ENDIF
              ELSE
                 ATIE(L) = 1.0
                 RISK(L) = REAL(NTOT - I)
                 L = L + 1
              ENDIF
           ELSE
              IC=IC+1
              ZC(IC)= XT(I)
           ENDIF
  300   CONTINUE
        LTOT = L - 1

        RETURN
        END

C
C      **********************************************************************
C      ********************** SUBROUTINE SORT1  *****************************
C      **********************************************************************
C
       SUBROUTINE SORT1(ID,X,Y,NTOT,INDEXX)
C
C      *       BUBBLE SORT PROGRAM                                          *
C      *       THIS PROGRAM ARRANGES OBSERVATIONS IN ASCENDING ORDER        *
C      *       ALSO IF THERE ARE TIED DATA POINTS, IT CHECKS THE CENSORING  *
C      *       STATUS AND ORDERS THEM SO THAT A DETECTED POINT COMES        *
C      *       FIRST.                                                       *
C      *                                                                    *
C      *       INPUT : INDEXX(I): POSITION INDICATOR                         *
C      *                 ID(I) : INDICATOR OF CENSORING                     *
C      *                 X(I): INDEPENDENT VARIABLE       *
C      *                 Y(I)  : DEPENDENT VARIABLE                         *
C      *                 NTOT  : NUMBER OF DATA POINTS                      *
C      *                                                                    *
C      *      OUTPUT :   ID, X, AND Y IN ASCENDING ORDER WITH INDEXX         *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       DIMENSION ID(NTOT),X(NTOT),Y(NTOT),INDEXX(NTOT)
C
C      *        SORTING IN Y, ASCENDING ORDER                               *
C
       DO 10 I=1,NTOT
          INDEXX(I)=I
   10  CONTINUE
C
       IF(NTOT.EQ.1) GOTO 100
C
       DO 90 I=1,NTOT
          J=NTOT-I+1
          JJ=J-1
          IF(JJ.GE.1) THEN   
C
             DO 80 K=1,JJ
                IF(Y(K).GT.Y(J)) THEN 

                   ID1=ID(J)
                   X1=X(J)
                   Y1=Y(J)
                   INDX=INDEXX(J)

                   ID(J)=ID(K)
                   X(J)=X(K)
                   Y(J)=Y(K)
                   INDEXX(J)=INDEXX(K)

                   ID(K)=ID1
                   X(K)=X1
                   Y(K)=Y1
                   INDEXX(K)=INDX
                ENDIF
  80         CONTINUE
          ENDIF
  90   CONTINUE
C
 100   RETURN
       END
