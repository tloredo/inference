C
C
C      *  ASURV:   SURVIVAL ANALYSIS PACKAGE FOR ASTRONOMERS                *
C      *                                                                    *
C      *  DEVELOPED BY:        TAKASHI ISOBE                                *
C      *                 CENTER FOR SPACE RESEARCH                          * 
C      *           THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY                *
C      *                                                                    *
C      *                      MICHAEL LAVALLEY                              *
C      *                  DEPARTMENT OF STATISTICS                          *
C      *               THE PENSYLVANIA STATE UNIVERSITY                     *
C      *       330A CLASSROOM BUILDING, UNIVERSITY PARK PA 16802            *
C      *                 INTERNET: MLV@STAT.PSU.EDU                         *
C      *                                                                    *
C      *                      ERIC FEIGELSON                                *
C      *             DEPARTMENT OF ASTRONOMY AND ASTROPHYSICS               * 
C      *                THE PENSYLVANIA STATE UNIVERSITY                    *
C      *              525 DAVEY LAB. UNIVERSITY PARK PA 16802               *
C      *                                                                    *
C      *  REV. 1.2  SECOND UPDATE   SUMMER 1992                             *
C      *                                                                    *
C      *     THIS PACKAGE IS WRITTEN TO PROVIDE SEVERAL                     *
C      *  SURVIVAL ANALYSIS METHODS WHICH ARE USEFUL IN ANALYZING           *
C      *  ASTRONOMICAL DATA. SURVIVAL ANALYSIS IS A GROUP OF STATISTICAL    *
C      *  METHODS WHICH TREAT PROBLEMS WITH CENSORED DATA (UPPER OR LOWER   *
C      *  LIMITS). THIS PACKAGE INCLUDES SOME TECHNIQUES DEVELOPED IN       *
C      *  IN OTHER FIELDS (E.G. ECONOMICS, ACTUARIAL SCIENCE, RELIABILITY   *
C      *  MATHEMATICS), AND A FEW METHODS DEVELOPED BY ASTRONOMERS.         *
C      *                                                                    *
C      *   THE METHODS PROVIDED IN THIS PACKAGE ARE :                       *
C      *                                                                    *
C      *   UNIVARIATE DISTRIBUTION :  KAPLAN-MEIER ESTIMATOR                *
C      *   TWO-SAMPLE TESTS        :  GEHAN TEST                            *
C      *                              LOGRANK TEST                          *
C      *                              PETO AND PETO TEST                    *
C      *                              PETO AND PRENTICE TEST                *
C      *   CORRELATION TESTS       :  COX PROPORTIONAL HAZARDS MODEL        *
C      *                              GENERALIZED KENDALL'S TAU (BHK METHOD)*
C      *                              GENERALIZED SPEARMAN'S RHO            *
C      *                                   (AKRITAS' METHOD)                *
C      *   LINEAR REGRESSIONS      :  EM ALGORITHM WITH NORMAL DISTRIBUTION *
C      *                              BUCKLEY-JAMES METHOD                  *
C      *                              TWO-DIMENSIONAL KAPLAN-MEIER          *
C      *                                  REGRESSION FOR DUAL-CENSORED DATA *
C      *                                                                    *
C      *                                                                    *
C      *   INPUTS                                                           *
C      *                                                                    *
C      *       IS0     :  IF 1 : UNIVARIATE PROBLEM                         *
C      *                     2 : CORRELATION/REGRESSION PROBLEM             *
C      *                     3 : EXIT                                       *
C      *                                                                    *
C      *   SUBROUTINES DATA1, UNIVAR, BIVAR                                 *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       CHARACTER*1 BLANK
C       OPEN(6,CARRIAGECONTROL='LIST',STATUS='OLD')
C
C
C
       PRINT *
       PRINT *,'    ***************************************************'
       PRINT *,'    *                                                 *'
       PRINT *,'    *                WELCOME TO ASURV                 *'
       PRINT *,'    *           SURVIVAL ANALYSIS PACKAGE             *' 
       PRINT *,'    *                FOR ASTRONOMERS                  *'
       PRINT *,'    *                                                 *'
       PRINT *,'    *                  DEVELOPED BY:                  *'
       PRINT *,'    *                  TAKASHI ISOBE                  *' 
       PRINT *,'    *         (CENTER FOR SPACE RESEARCH, MIT)        *'
       PRINT *,'    *                 MICHAEL LAVALLEY                *'
       PRINT *,'    *         (DEPT. OF STATISTICS, PENN STATE)       *'
       PRINT *,'    *                  ERIC FEIGELSON                 *'
       PRINT *,'    * (DEPT. OF ASTRONOMY & ASTROPHYSICS, PENN STATE) *'
       PRINT *,'    *                                                 *'
       PRINT *,'    *                                                 *'
       PRINT *,'    *              REV 1.2  SUMMER 1992               *'
       PRINT *,'    ***************************************************'
       PRINT *
       PRINT *
       PRINT *
       PRINT *
       PRINT *,'              (CARRIAGE RETURN TO CONTINUE) '
       READ(5,50) BLANK
   50  FORMAT(A1)
       PRINT *
C
C      *           START CONVERSATION WITH THE USER                         *
C
       PRINT *
       PRINT *
       PRINT *
  100  PRINT *,'                          MENU  '
       PRINT *
       PRINT *
       PRINT *,'       UNIVARIATE DATA           BIVARIATE DATA '
       PRINT *
       PRINT *
       PRINT *,'     DISTRIBUTION FUNCTION       CORRELATION '
       PRINT *,'   1 KAPLAN-MEIER ESTIMATOR    1 COX REGRESSION '
       PRINT *,'                               2 GEN. KENDALL TAU'
       PRINT *,'                               3 GEN. SPEARMAN RHO'
       PRINT *
       PRINT *
       PRINT *,'     TWO-SAMPLE TESTS            LINEAR REGRESSION '
       PRINT *,'   1 GEHAN TESTS               1 EM ALGORITHM WITH  '
       PRINT *,'   2 LOGRANK TEST                 GAUSSIAN RESIDUALS ' 
       PRINT *,'   3 PETO AND PETO TEST        2 BUCKLEY-JAMES METHOD ' 
       PRINT *,'   4 PETO AND PRENTICE TEST       WITH KM RESIDUALS '
       PRINT *,'                               3 SCHMITT METHOD FOR '
       PRINT *,'                                  DUAL CENSORED DATA ' 
       PRINT *
       PRINT *
       PRINT *
       PRINT *,'            (CARRIAGE RETURN TO CONTINUE) '
       READ(5,50) BLANK
C
       PRINT *
C
C      *  CHOICE : UNIVARIATE PROBLEM OR CORRELATION/REGRESSION PROBLEM     *
C
       PRINT *
       PRINT *,'    SELECT DATA TYPE: ' 
       PRINT *,'     1 UNIVARIATE DATA '
       PRINT *,'     2 BIVARIATE DATA ' 
       PRINT *,'     3 EXIT '
  200  WRITE(6,210)
  210  FORMAT(' CHOICE ? ')
C  210  FORMAT('          CHOICE ? ',$)
C
       CALL DATA1(IS0)
C
       IF((IS0.EQ.1).OR.(IS0.EQ.2).OR.(IS0.EQ.3)) GOTO 300
       PRINT *,'PLEASE TYPE ONCE MORE'
       GOTO 200
C
  300  IBACK=0
       IF(IS0.EQ.1) CALL UNIVAR(IBACK)
       IF(IS0.EQ.2) CALL BIVAR(IBACK)
       IF(IS0.EQ.3) STOP
C
       IF(IBACK.EQ.1) GOTO 100
       STOP
       END

