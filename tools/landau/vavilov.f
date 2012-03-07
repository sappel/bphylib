      FUNCTION DISVAV (X,I)      ! mod R.W. Hasse 25.4.2005

                                                                        AAUJ0055
C     DISVAV(X,I) COMPUTES FOR I=0 THE VALUE OF THE DENSITY, FOR I=1    AAUJ0056
C     THE VALUE OF THE CUMULATIVE DISTRIBUTION FUNCTION OF THE          AAUJ0057
C     VAVILOV DISTRIBUTION AT THE POINT X                               AAUJ0058
                                                                        AAUJ0059
      COMMON /VAVILI/ T0,T1,T,OMEGA                                     AAUJ0060
      COMMON /VAVILO/ A(155),B(155),N                                   AAUJ0061
                                                                        AAUJ0062
      IF(X .LT. T0) GO TO 3                                             AAUJ0063
      IF(X .GT. T1) GO TO 4                                             AAUJ0064
                                                                        AAUJ0065
      Y=X-T0                                                            AAUJ0066
                                                                        AAUJ0067
      TRR=0.                                                            AAUJ0068
      IF(I .EQ. 1) TRR=Y/T                                              AAUJ0069
                                                                        AAUJ0070
      Z=OMEGA*Y-3.1415926535898                                         AAUJ0071
      COF=2.0*COS(Z)                                                    AAUJ0072
      A1=0.                                                             AAUJ0073
      A0=A(1)                                                           AAUJ0074
      NN=N-1                                                            AAUJ0075
      DO 1 K = 2,N                                                      AAUJ0076
      A2=A1                                                             AAUJ0077
      A1=A0                                                             AAUJ0078
    1 A0=A(K)+COF*A1-A2                                                 AAUJ0079
                                                                        AAUJ0080
      TRR=TRR+0.5*(A0-A2)                                               AAUJ0081
      A1=0.                                                             AAUJ0082
      A0=B(1)                                                           AAUJ0083
      DO 2 K = 2,NN                                                     AAUJ0084
      A2=A1                                                             AAUJ0085
      A1=A0                                                             AAUJ0086
    2 A0=B(K)+COF*A1-A2                                                 AAUJ0087
                                                                        AAUJ0088
      DISVAV=TRR+A0*SIN(Z)                                              AAUJ0089
      RETURN                                                            AAUJ0090
                                                                        AAUJ0091
    3 DISVAV=0.0                                                        AAUJ0092
      RETURN                                                            AAUJ0093
                                                                        AAUJ0094
    4 DISVAV=I                                                          AAUJ0095
      RETURN                                                            AAUJ0096
                                                                        AAUJ0097
      END                                                               AAUJ0098

c-----------------------------------------------------------------------

      FUNCTION DINVAV(X)                                                AAUJ0099
                                                                        AAUJ0100
C     DINVAV(X) COMPUTES THE VALUE OF THE INVERSE OF THE CONDITIONAL    AAUJ0101
C     CUMULATIVE DISTRIBUTION FUNCTION OF THE VAVILOV DISTRIBUTION      AAUJ0102
C     AT THE POINT X                                                    AAUJ0103
                                                                        AAUJ0104
      COMMON /VAVILA/ B(200),S,TT                                       AAUJ0105
      COMMON /VAVILI/ T0,T1,T,OMEGA                                     AAUJ0106
                                                                        AAUJ0107
      Z=3.1415926535898*X                                               AAUJ0108
      COF=2.*COS(Z)                                                     AAUJ0109
      C1=0.                                                             AAUJ0110
      C0=B(1)                                                           AAUJ0111
                                                                        AAUJ0112
      DO 1 K = 2,200                                                    AAUJ0113
      C2=C1                                                             AAUJ0114
      C1=C0                                                             AAUJ0115
    1 C0=B(K)+COF*C1-C2                                                 AAUJ0116
                                                                        AAUJ0117
      DINVAV=T0+TT*X+C0*SIN(Z)                                          AAUJ0118
      RETURN                                                            AAUJ0119
                                                                        AAUJ0120
      END

c-----------------------------------------------------------------------

      SUBROUTINE COEDIS(RKA,BE2,I,J)                                    AAUJ0122
                                                                        AAUJ0123
C     COEDIS COMPUTES THE ENDPOINTS T- AND T+ OF THE SUPPORT OF         AAUJ0124
C     DISVAV(X,0).IT ALSO COMPUTES THE FOURIER COEFFICIENTS OF DISVAV(X,AAUJ0125
                                                                        AAUJ0126
      COMMON /VAVILI/ T0,T1,T,OMEGA                                     AAUJ0127
      COMMON /VAVILO/ A(155),B(155),N                                   AAUJ0128
      COMMON /FORFCN/ SS,LFCN                                           AAUJ0129
      DIMENSION XP(8),XQ(6)                                             AAUJ0130
      DATA E,PI,RG /5E-4, 3.1415926535898, 0.5772156649015/             AAUJ0131
      DATA (XP(I),I=1,8)                                                AAUJ0132
     1 /9.29, 2.47, 0.89, 0.36, 0.15, 0.07, 0.03, 0.02/                 AAUJ0133
      DATA (XQ(I),I=1,6)                                                AAUJ0134
     1 /0.012, 0.03, 0.08, 0.26, 0.87, 3.83/                            AAUJ0135
                                                                        AAUJ0136
      LU=IABS(J)                                                        AAUJ0137
      IF (RKA .LT. 0.01 .OR. RKA .GT. 10.0) GOTO 6                      AAUJ0138
      IF (BE2 .LT. 0.0  .OR. BE2 .GT.  1.0) GOTO 8                      AAUJ0139
                                                                        AAUJ0140
      Z=1.-BE2*(1.-RG)-ALOG(E)/RKA                                      AAUJ0141
      T0=(ALOG(E)/RKA-(1.+BE2*RG)-Z*ALOG(RKA)-(Z+BE2)*(ALOG(Z)          AAUJ0142
     1  +EXPINT(Z))+EXP(-Z))/Z                                          AAUJ0143
      DO 1 L = 1,8                                                      AAUJ0144
      IF(RKA .GE. XP(L)) GO TO 11                                       AAUJ0145
    1 CONTINUE                                                          AAUJ0146
      L=9                                                               AAUJ0147
   11 P=-FLOAT(L)-0.5                                                   AAUJ0148
      DO 2 L = 1,6                                                      AAUJ0149
      IF(RKA .LE. XQ(L)) GO TO 22                                       AAUJ0150
    2 CONTINUE                                                          AAUJ0151
      L=7                                                               AAUJ0152
   22 Q=FLOAT(L)-7.5                                                    AAUJ0153
      LFCN=3                                                            AAUJ0154
      CALL RZERO(P,Q,X,RKA,BE2,LU)                                      AAUJ0155
      T1=(ALOG(E)/RKA-(1.+BE2*RG))/X-ALOG(RKA)-(1.+BE2/X)*(ALOG(ABS(X)) AAUJ0156
     1  +EXPINT(X))+EXP(-X)/X                                           AAUJ0157
                                                                        AAUJ0158
      IF(J .GT. 0) Continue ! WRITE(J,10) T0,T1
      T=T1-T0                                                           AAUJ0160
      OMEGA=2.0*PI/T                                                    AAUJ0161
      LFCN=1                                                            AAUJ0162
      CALL RZERO(5.,155.,X,RKA,BE2,LU)                                  AAUJ0163
      N=X+1.                                                            AAUJ0164
                                                                        AAUJ0165
      D=EXP(RKA*(1.+BE2*(RG-ALOG(RKA))))                                AAUJ0166
      A(N)=0.                                                           AAUJ0167
      IF(I .EQ. 0) A(N)=OMEGA/PI                                        AAUJ0168
      N1=N-1                                                            AAUJ0169
      Q=-1.                                                             AAUJ0170
                                                                        AAUJ0171
      DO 3 K = 1,N1                                                     AAUJ0172
      L=N-K                                                             AAUJ0173
      X=OMEGA*FLOAT(K)                                                  AAUJ0174
      X1=X/RKA                                                          AAUJ0175
      C1=ALOG(X)-CSINTL(X1,2)                                           AAUJ0176
      C2=CSINTL(X1,1)                                                   AAUJ0177
      C3=SIN(X1)                                                        AAUJ0178
      C4=COS(X1)                                                        AAUJ0179
      F1=RKA*(BE2*C1-C4)-X*C2                                           AAUJ0180
      F2=X*C1+RKA*(C3+BE2*C2)+T0*X                                      AAUJ0181
      D1=Q*D*EXP(F1)/PI                                                 AAUJ0182
      HS=D1*SIN(F2)                                                     AAUJ0183
      HC=D1*COS(F2)                                                     AAUJ0184
      IF(I .EQ. 0) GO TO 4                                              AAUJ0185
      A(L)=HS/FLOAT(K)                                                  AAUJ0186
      B(L)=HC/FLOAT(K)                                                  AAUJ0187
      A(N)=A(N)-2.0*Q*A(L)                                              AAUJ0188
      GO TO 3                                                           AAUJ0189
    4 A(L)=HC*OMEGA                                                     AAUJ0190
      B(L)=-HS*OMEGA                                                    AAUJ0191
    3 Q=-Q                                                              AAUJ0192
      RETURN                                                            AAUJ0193
                                                                        AAUJ0194
    6 Continue ! WRITE(LU,7) RKA
      RETURN                                                            AAUJ0196
                                                                        AAUJ0197
    8 Return ! WRITE(LU,9) BE2
    7 FORMAT(/10X,7HKAPPA =,E10.3,15H - OUT OF RANGE)                   AAUJ0199
    9 FORMAT(/10X,9HBETA**2 =,E10.3,15H - OUT OF RANGE)                 AAUJ0200
   10 FORMAT(10X,4HT- =,F8.3,10X,4HT+ =,F8.3)                           AAUJ0201
      RETURN                                                            AAUJ0202
                                                                        AAUJ0203
      END

c-----------------------------------------------------------------------

      SUBROUTINE COEDIN(RKA,BE2,J)                                      AAUJ0205
                                                                        AAUJ0206
C     COEDIN COMPUTES THE FOURIER COEFFICIENTS FOR THE INVERSE OF       AAUJ0207
C     THE CONDITIONAL CUMULATIVE DISTRIBUTION FUNCTION OF THE           AAUJ0208
C     VAVILOV DISTRIBUTION                                              AAUJ0209
                                                                        AAUJ0210
      DIMENSION C(1001)                                                 AAUJ0211
      COMMON /VAVILA/ B(200),S,TT                                       AAUJ0212
      COMMON /VAVILI/ T0,T1,T,OMEGA                                     AAUJ0213
      COMMON /FORFCN/ SS,LFCN                                           AAUJ0214
      DATA PI/3.1415926535898/                                          AAUJ0215
                                                                        AAUJ0216
      LU=IABS(J)                                                        AAUJ0217
      CALL COEDIS(RKA,BE2,1,J)                                          AAUJ0218
                                                                        AAUJ0219
      SS=0.99                                                           AAUJ0220
      IF(RKA .GE. 0.04) SS=0.995                                        AAUJ0221
      LFCN=2                                                            AAUJ0222
      CALL RZERO(T0,T1,TR,RKA,BE2,LU)                                   AAUJ0223
      IF(J .GT. 0) Continue ! WRITE(J,4) TR
      S=DISVAV(TR,1)                                                    AAUJ0225
      TT=TR-T0                                                          AAUJ0226
      STEP=TT/1000.                                                     AAUJ0227
      STEP2=2.0*STEP                                                    AAUJ0228
                                                                        AAUJ0229
      DO 1 I = 1,1001                                                   AAUJ0230
    1 C(I)=DISVAV(T0+STEP*FLOAT(I-1),1)                                 AAUJ0231
                                                                        AAUJ0232
      DO 2 K = 1,200                                                    AAUJ0233
      K1=201-K                                                          AAUJ0234
      Z1=PI*FLOAT(K)                                                    AAUJ0235
      Z=Z1/S                                                            AAUJ0236
      BK=0.5*(COS(Z*C(1))+COS(Z*C(1001)))                               AAUJ0237
      DO 3 I = 2,1000                                                   AAUJ0238
    3 BK=BK+COS(Z*C(I))                                                 AAUJ0239
    2 B(K1)=STEP2*BK/Z1                                                 AAUJ0240
    4 FORMAT(10X,9HT PRIME =,F10.4)                                     AAUJ0241
      RETURN                                                            AAUJ0242
                                                                        AAUJ0243
      END                                                               AAUJ0244

c-----------------------------------------------------------------------

      FUNCTION FCN(X,RKA,BE2)                                           AAUJ0246
                                                                        AAUJ0247
      COMMON /VAVILI/ T0,T1,T,OMEGA                                     AAUJ0248
      COMMON /FORFCN/ SS,LFCN                                           AAUJ0249
      DATA E,PI,RG /5E-4, 3.1415926535898, 0.5772156649015/             AAUJ0250
                                                                        AAUJ0251
      GO TO (1,2,3), LFCN                                               AAUJ0252
                                                                        AAUJ0253
C     FOR LFCN=1 FCN IS USED TO DETERMINE THE NUMBER N OF FOURIER       AAUJ0254
C     COEFFICIENTS IN DISVAV(X,I)                                       AAUJ0255
                                                                        AAUJ0256
    1 RN=5.                                                             AAUJ0257
      IF(RKA .GE. 0.07) RN=8.                                           AAUJ0258
      FCN=RKA*(2.+BE2*RG)-1.596312592                                   AAUJ0259
     1   +BE2*RKA*ALOG(OMEGA*X/RKA)-0.5*OMEGA*PI*X+2.302585093*RN       AAUJ0260
      RETURN                                                            AAUJ0261
                                                                        AAUJ0262
C     FOR LFCN=2 FCN IS USED TO DETERMINE THE SS-QUANTILE OF THE        AAUJ0263
C     VAVILOV DISTRIBUTION                                              AAUJ0264
                                                                        AAUJ0265
    2 FCN=DISVAV(X,1)-SS                                                AAUJ0266
      RETURN                                                            AAUJ0267
                                                                        AAUJ0268
C     FOR LFCN=3 FCN IS USED TO DETERMINE THE RIGHT-HAND ENDPOINT T+    AAUJ0269
C     OF THE SUPPORT OF DISVAV(X,0)                                     AAUJ0270
                                                                        AAUJ0271
    3 FCN=1.-BE2*(1.-RG)-ALOG(E)/RKA-X+BE2*(ALOG(ABS(X))+EXPINT(X))     AAUJ0272
     1    -(1.-BE2)*EXP(-X)                                             AAUJ0273
      RETURN                                                            AAUJ0274
                                                                        AAUJ0275
      END

c-----------------------------------------------------------------------

      FUNCTION EXPINT(X)                                                AAUJ0277
                                                                        AAUJ0278
C     EXPINT(X) COMPUTES THE VALUE OF THE EXPONENTIAL INTEGRAL          AAUJ0279
C     FOR ANY VALUE OF X                                                AAUJ0280
                                                                        AAUJ0281
      DIMENSION P1(5),Q1(4),P2(7),Q2(7),P3(6),Q3(6),P4(8),Q4(8)         AAUJ0282
      DIMENSION A1(8),B1(7),A2(8),B2(7),A3(6),B3(5)                     AAUJ0283
                                                                        AAUJ0284
      DATA (P1(I),I=1,5)                                                AAUJ0285
     1/-4.34981 43832 95212E+2, +4.25696 82638 59170E+2,                AAUJ0286
     2 +2.92525 18866 92055E+2, +3.98941 53870 32107E+1,                AAUJ0287
     3 +4.29312 52343 20973E+0/                                         AAUJ0288
      DATA (Q1(I),I=1,4)                                                AAUJ0289
     1/+7.53585 64359 84293E+2, +5.68052 52718 98696E+2,                AAUJ0290
     2 +1.50950 38744 25131E+2, +1.88992 88395 00297E+1/                AAUJ0291
      DATA (P2(I),I=1,7)                                                AAUJ0292
     1/+4.65627 10797 50957E-7, +9.99979 57705 15950E-1,                AAUJ0293
     2 +9.04161 55694 63287E+0, +2.43784 08879 13167E+1,                AAUJ0294
     3 +2.30192 55939 13335E+1, +6.90522 52278 44436E+0,                AAUJ0295
     4 +4.30967 83946 93888E-1/                                         AAUJ0296
      DATA (Q2(I),I=1,7)                                                AAUJ0297
     1/+1.00000 00000 00000E+0, +1.00411 64382 90545E+1,                AAUJ0298
     2 +3.24264 21069 51381E+1, +4.12807 84189 14243E+1,                AAUJ0299
     3 +2.04494 78501 37942E+1, +3.31909 21359 33016E+0,                AAUJ0300
     4 +1.03400 13040 48740E-1/                                         AAUJ0301
      DATA (P3(I),I=1,6)                                                AAUJ0302
     1/-9.99999 99990 36009E-1, -1.96304 08535 93865E+1,                AAUJ0303
     2 -1.19557 61038 37175E+2, -2.54376 33976 88951E+2,                AAUJ0304
     3 -1.47982 19500 50448E+2, -2.39099 64453 13573E+0/                AAUJ0305
      DATA (Q3(I),I=1,6)                                                AAUJ0306
     1/+1.00000 00000 00000E+0, +2.16304 08494 23774E+1,                AAUJ0307
     2 +1.56818 43364 53856E+2, +4.62230 27156 14783E+2,                AAUJ0308
     3 +5.30685 09610 81167E+2, +1.77600 70940 35063E+2/                AAUJ0309
      DATA (P4(I),I=1,8)                                                AAUJ0310
     1/-8.66937 33995 10696E+0, -5.49142 26552 10851E+2,                AAUJ0311
     2 -4.21001 61535 70699E+3, -2.49301 39345 86476E+5,                AAUJ0312
     3 -1.19623 66934 92469E+5, -2.21744 62775 88454E+7,                AAUJ0313
     4 +3.89280 42131 12014E+6, -3.91546 07380 90955E+8/                AAUJ0314
      DATA (Q4(I),I=1,8)                                                AAUJ0315
     1/+3.41718 75000 00000E+1, -1.60708 92658 72209E+3,                AAUJ0316
     2 +3.57300 29805 85081E+4, -4.83547 43616 21635E+5,                AAUJ0317
     3 +4.28559 62461 17490E+6, -2.49033 37574 05403E+7,                AAUJ0318
     4 +8.91925 76757 56121E+7, -1.65254 29972 52109E+8/                AAUJ0319
      DATA (A1(I),I=1,8)                                                AAUJ0320
     1/+1.00443 10922 80779E+0, -4.32531 13287 81346E+1,                AAUJ0321
     2 +6.01217 99083 00805E+1, -3.31842 53199 72211E+1,                AAUJ0322
     3 +2.50762 81129 35598E+1, +9.30816 38566 21651E+0,                AAUJ0323
     4 -2.19010 23385 48809E+1, -2.18086 38152 07237E+0/                AAUJ0324
      DATA (B1(I),I=1,7)                                                AAUJ0325
     1/+5.27468 85196 29079E-1, +2.73624 11988 93281E+3,                AAUJ0326
     2 +1.43256 73812 19376E+1, +1.00367 43951 67258E+3,                AAUJ0327
     3 -6.25041 16167 18755E+0, +3.00892 64837 29152E+2,                AAUJ0328
     4 +3.93707 70185 27150E+0/                                         AAUJ0329
      DATA (A2(I),I=1,8)                                                AAUJ0330
     1/+9.99994 29607 47083E-1, -1.95022 32128 96598E+0,                AAUJ0331
     2 +1.75656 31546 96144E+0, +1.79601 68876 92516E+1,                AAUJ0332
     3 -3.23467 33030 54035E+1, -8.28561 99414 06413E+0,                AAUJ0333
     4 -1.86545 45488 33988E+1, -3.48334 65360 28526E+0/                AAUJ0334
      DATA (B2(I),I=1,7)                                                AAUJ0335
     1/+1.00083 86740 26391E+0, -3.43942 26689 98700E+0,                AAUJ0336
     2 +2.89516 72792 51351E+1, +7.60761 14800 77346E+2,                AAUJ0337
     3 +2.57776 38423 84399E+1, +5.72837 19383 73237E+1,                AAUJ0338
     4 +6.95000 65588 74340E+1/                                         AAUJ0339
      DATA (A3(I),I=1,6)                                                AAUJ0340
     1/+1.00000 00000 70443E+0, -3.00000 77799 35772E+0,                AAUJ0341
     2 -5.02233 17461 85109E+0, -9.14830 08216 73641E+0,                AAUJ0342
     3 -1.01047 90815 76032E+1, -2.77809 28934 43810E+1/                AAUJ0343
      DATA (B3(I),I=1,5)                                                AAUJ0344
     1/+1.99999 99428 26009E+0, -2.99901 18065 26193E+0,                AAUJ0345
     2 -7.18975 18395 04450E+0, +2.72761 00778 77917E+0,                AAUJ0346
     3 +1.22399 93926 82269E+2/                                         AAUJ0347
      DATA X0 /0.37250 74107 81367/                                     AAUJ0348
                                                                        AAUJ0349
      IF(X .GT. 4.0) GO TO 1                                            AAUJ0350
      IF(X .GT. 1.0) GO TO 2                                            AAUJ0351
      IF(X .GT. 0.0) GO TO 3                                            AAUJ0352
      IF(X .EQ. 0.0) GO TO 4                                            AAUJ0353
      Y=-X                                                              AAUJ0354
      IF(X .GT. -6.0) GO TO 5                                           AAUJ0355
      V=EXP(Y)/X                                                        AAUJ0356
      IF(X .GT. -24.0) GO TO 7                                          AAUJ0357
                                                                        AAUJ0358
      EXPINT=V*(1.0+(A3(1)+B3(1)/(A3(2)+Y+B3(2)/(A3(3)+Y+B3(3)/(A3(4)+Y+AAUJ0359
     1 B3(4)/(A3(5)+Y+B3(5)/(A3(6)+Y))))))/Y)                           AAUJ0360
      RETURN                                                            AAUJ0361
                                                                        AAUJ0362
    7 EXPINT=V*(A2(1)+B2(1)/(A2(2)+Y+B2(2)/(A2(3)+Y+B2(3)/(A2(4)+Y+B2(4)AAUJ0363
     1 /(A2(5)+Y+B2(5)/(A2(6)+Y+B2(6)/(A2(7)+Y+B2(7)/(A2(8)+Y))))))))   AAUJ0364
      RETURN                                                            AAUJ0365
                                                                        AAUJ0366
    6 EXPINT=V*(A1(1)+B1(1)/(A1(2)+Y+B1(2)/(A1(3)+Y+B1(3)/(A1(4)+Y+B1(4)AAUJ0367
     1 /(A1(5)+Y+B1(5)/(A1(6)+Y+B1(6)/(A1(7)+Y+B1(7)/(A1(8)+Y))))))))   AAUJ0368
      RETURN                                                            AAUJ0369
                                                                        AAUJ0370
    5 V=0.666666666666667*Y-2.0                                         AAUJ0371
      BP=0.                                                             AAUJ0372
      BQ=0.                                                             AAUJ0373
      DP=P4(1)                                                          AAUJ0374
      DQ=Q4(1)                                                          AAUJ0375
      DO 15 I = 2,8                                                     AAUJ0376
      AP=BP                                                             AAUJ0377
      BP=DP                                                             AAUJ0378
      AQ=BQ                                                             AAUJ0379
      BQ=DQ                                                             AAUJ0380
      DP=P4(I)-AP+V*BP                                                  AAUJ0381
   15 DQ=Q4(I)-AQ+V*BQ                                                  AAUJ0382
      EXPINT=-ALOG(Y/X0)-(Y-X0)*(DP-AP)/(DQ-AQ)                         AAUJ0383
      RETURN                                                            AAUJ0384
                                                                        AAUJ0385
    4 EXPINT=0.                                                         AAUJ0386
      RETURN                                                            AAUJ0387
                                                                        AAUJ0388
    3 EXPINT=-ALOG(X)+                                                  AAUJ0389
     1       (P1(1)+X*(P1(2)+X*(P1(3)+X*(P1(4)+X*P1(5)))))/             AAUJ0390
     2       (Q1(1)+X*(Q1(2)+X*(Q1(3)+X*(Q1(4)+X))))                    AAUJ0391
      RETURN                                                            AAUJ0392
                                                                        AAUJ0393
    2 Y=1.0/X                                                           AAUJ0394
      EXPINT=EXP(-X)*                                                   AAUJ0395
     1(P2(1)+Y*(P2(2)+Y*(P2(3)+Y*(P2(4)+Y*(P2(5)+Y*(P2(6)+Y*P2(7)))))))/AAUJ0396
     2(Q2(1)+Y*(Q2(2)+Y*(Q2(3)+Y*(Q2(4)+Y*(Q2(5)+Y*(Q2(6)+Y*Q2(7))))))) AAUJ0397
      RETURN                                                            AAUJ0398
                                                                        AAUJ0399
    1 Y=1.0/X                                                           AAUJ0400
      EXPINT=EXP(-X)*Y*(1.0+Y*                                          AAUJ0401
     1       (P3(1)+Y*(P3(2)+Y*(P3(3)+Y*(P3(4)+Y*(P3(5)+Y*P3(6))))))/   AAUJ0402
     2       (Q3(1)+Y*(Q3(2)+Y*(Q3(3)+Y*(Q3(4)+Y*(Q3(5)+Y*Q3(6)))))))   AAUJ0403
      RETURN                                                            AAUJ0404
                                                                        AAUJ0405
      END                                                               AAUJ0406

c-----------------------------------------------------------------------

      FUNCTION CSINTL(X,L)                                              AAUJ0407
                                                                        AAUJ0408
C     CSINTL(X,L) COMPUTES FOR L=1 THE SINE INTEGRAL AND FOR            AAUJ0409
C     L=2 THE COSINE INTEGRAL FOR ANY VALUE OF X                        AAUJ0410
                                                                        AAUJ0411
      GO TO (11,12), L                                                  AAUJ0412
                                                                        AAUJ0413
   11 IF(ABS(X) .GE. 16.0) GO TO 1                                      AAUJ0414
      Z=X/16.0                                                          AAUJ0415
      Y=4.0*Z**2-2.0                                                    AAUJ0416
      B=     +0.00000 00000 00007                                       AAUJ0417
      A=Y*B  -0.00000 00000 00185                                       AAUJ0418
      B=Y*A-B+0.00000 00000 04185                                       AAUJ0419
      A=Y*B-A-0.00000 00000 84710                                       AAUJ0420
      B=Y*A-B+0.00000 00015 22370                                       AAUJ0421
      A=Y*B-A-0.00000 00241 00076                                       AAUJ0422
      B=Y*A-B+0.00000 03329 88589                                       AAUJ0423
      A=Y*B-A-0.00000 39729 08746                                       AAUJ0424
      B=Y*A-B+0.00004 04202 71419                                       AAUJ0425
      A=Y*B-A-0.00034 54691 55569                                       AAUJ0426
      B=Y*A-B+0.00243 62214 04749                                       AAUJ0427
      A=Y*B-A-0.01386 74455 89417                                       AAUJ0428
      B=Y*A-B+0.06203 36794 32003                                       AAUJ0429
      A=Y*B-A-0.21126 37809 76555                                       AAUJ0430
      B=Y*A-B+0.53014 88479 16522                                       AAUJ0431
      A=Y*B-A-0.96832 22369 87086                                       AAUJ0432
      B=Y*A-B+1.38930 87711 71888                                       AAUJ0433
      A=Y*B-A-1.92656 50911 50656                                       AAUJ0434
      B=Y*A-B+2.77875 63817 42663                                       AAUJ0435
      A=Y*B-A-4.06398 08449 11986                                       AAUJ0436
      A=Y*A-B+8.10585 29553 61245                                       AAUJ0437
      CSINTL=Z*0.5*(A-B)                                                AAUJ0438
      RETURN                                                            AAUJ0439
                                                                        AAUJ0440
    1 Z=16.0/X                                                          AAUJ0441
      G=Z**2                                                            AAUJ0442
      Y=4.0*G-2.0                                                       AAUJ0443
      B=     +0.00000 00000 00002                                       AAUJ0444
      A=Y*B  -0.00000 00000 00014                                       AAUJ0445
      B=Y*A-B+0.00000 00000 00107                                       AAUJ0446
      A=Y*B-A-0.00000 00000 00964                                       AAUJ0447
      B=Y*A-B+0.00000 00000 10308                                       AAUJ0448
      A=Y*B-A-0.00000 00001 36096                                       AAUJ0449
      B=Y*A-B+0.00000 00023 56196                                       AAUJ0450
      A=Y*B-A-0.00000 00586 70317                                       AAUJ0451
      B=Y*A-B+0.00000 24537 55677                                       AAUJ0452
      A=Y*B-A-0.00023 37560 41393                                       AAUJ0453
      A=Y*A-B+0.12452 74580 57854                                       AAUJ0454
      F=Z*0.5*(A-B)                                                     AAUJ0455
      B=     +0.00000 00000 00002                                       AAUJ0456
      A=Y*B  -0.00000 00000 00012                                       AAUJ0457
      B=Y*A-B+0.00000 00000 00087                                       AAUJ0458
      A=Y*B-A-0.00000 00000 00717                                       AAUJ0459
      B=Y*A-B+0.00000 00000 06875                                       AAUJ0460
      A=Y*B-A-0.00000 00000 79604                                       AAUJ0461
      B=Y*A-B+0.00000 00011 69202                                       AAUJ0462
      A=Y*B-A-0.00000 00234 68225                                       AAUJ0463
      B=Y*A-B+0.00000 07249 95950                                       AAUJ0464
      A=Y*B-A-0.00004 26441 82622                                       AAUJ0465
      A=Y*A-B+0.00772 57121 93407                                       AAUJ0466
      G=G*0.5*(A-B)                                                     AAUJ0467
                                                                        AAUJ0468
      GO TO (2,3), L                                                    AAUJ0469
                                                                        AAUJ0470
    2 B=1.57079 63267 94897                                             AAUJ0471
      IF(X .LT. 0.0) B=-B                                               AAUJ0472
      CSINTL=B-F*COS(X)-G*COS(X)                                        AAUJ0473
      RETURN                                                            AAUJ0474
                                                                        AAUJ0475
    3 CSINTL=F*SIN(X)-G*COS(X)                                          AAUJ0476
      RETURN                                                            AAUJ0477
                                                                        AAUJ0478
   12 F=ABS(X)                                                          AAUJ0479
      IF(F .GE. 16.0) GO TO 1                                           AAUJ0480
      Z=X/16.0                                                          AAUJ0481
      G=Z**2                                                            AAUJ0482
      Y=4.0*G-2.0                                                       AAUJ0483
      B=     + 0.00000 00000 00003                                      AAUJ0484
      A=Y*B  - 0.00000 00000 00078                                      AAUJ0485
      B=Y*A-B+ 0.00000 00000 01861                                      AAUJ0486
      A=Y*B-A- 0.00000 00000 40066                                      AAUJ0487
      B=Y*A-B+ 0.00000 00007 69379                                      AAUJ0488
      A=Y*B-A- 0.00000 00130 84020                                      AAUJ0489
      B=Y*A-B+ 0.00000 01954 37144                                      AAUJ0490
      A=Y*B-A- 0.00000 25401 25611                                      AAUJ0491
      B=Y*A-B+ 0.00002 84169 45498                                      AAUJ0492
      A=Y*B-A- 0.00027 02123 31184                                      AAUJ0493
      B=Y*A-B+ 0.00215 20467 52074                                      AAUJ0494
      A=Y*B-A- 0.01411 08652 53535                                      AAUJ0495
      B=Y*A-B+ 0.07467 82552 94576                                      AAUJ0496
      A=Y*B-A- 0.31208 09248 25428                                      AAUJ0497
      B=Y*A-B+ 1.00866 07873 58110                                      AAUJ0498
      A=Y*B-A- 2.49750 50885 39025                                      AAUJ0499
      B=Y*A-B+ 4.86202 23485 00627                                      AAUJ0500
      A=Y*B-A- 8.10790 39705 62531                                      AAUJ0501
      B=Y*A-B+12.74187 08697 58071                                      AAUJ0502
      A=Y*B-A-19.38612 40966 07770                                      AAUJ0503
      A=Y*A-B+29.98517 87356 26818                                      AAUJ0504
      CSINTL=ALOG(F)-G*0.5*(A-B)+0.57721 56649 01533                    AAUJ0505
      RETURN                                                            AAUJ0506
                                                                        AAUJ0507
      END                                                               AAUJ0508

c-----------------------------------------------------------------------

      SUBROUTINE RZERO(A,B,X,RKA,BE2,LU)                                AAUJ0509
                                                                        AAUJ0510
C     RZERO SEARCHES FOR THE ROOT OF THE EQUATION FCN=0 IN THE INTERVAL AAUJ0511
C     (A,B)                                                             AAUJ0512
      COMMON /FORFCN/ SS,LFCN                                           AAUJ0513
      DATA E,EPSI,MAXFUN/1E-9,1E-5,100/                                 AAUJ0514
      MC=0                                                              AAUJ0515
                                                                        AAUJ0516
      XA=AMIN1(A,B)                                                     AAUJ0517
      XB=AMAX1(A,B)                                                     AAUJ0518
      FA=FCN(A,RKA,BE2)                                                 AAUJ0519
      MC=MC+1                                                           AAUJ0520
      FB=FCN(B,RKA,BE2)                                                 AAUJ0521
      IF(FA*FB .GT. 0.0) GO TO 16                                       AAUJ0522
      MC=MC+1                                                           AAUJ0523
                                                                        AAUJ0524
    4 X=0.5*(XA+XB)                                                     AAUJ0525
      R=X-XA                                                            AAUJ0526
      EE=ABS(X)+E                                                       AAUJ0527
      IF(R .LE. EE*EPSI) GO TO 18                                       AAUJ0528
      F1=FA                                                             AAUJ0529
      X1=XA                                                             AAUJ0530
      F2=FB                                                             AAUJ0531
      X2=XB                                                             AAUJ0532
    1 MC=MC+1                                                           AAUJ0533
      G=FCN(X,RKA,BE2)                                                  AAUJ0534
      IF(MC .GT. MAXFUN) GO TO 17                                       AAUJ0535
      FX=G                                                              AAUJ0536
                                                                        AAUJ0537
      IF(FX*FA .GT. 0.0) GO TO 2                                        AAUJ0538
      FB=FX                                                             AAUJ0539
      XB=X                                                              AAUJ0540
      GO TO 3                                                           AAUJ0541
    2 XA=X                                                              AAUJ0542
      FA=FX                                                             AAUJ0543
                                                                        AAUJ0544
C     PARABOLA ITERATION                                                AAUJ0545
                                                                        AAUJ0546
    3 IF((X1-X2)*(X2-X)*(X1-X) .EQ. 0.0) GO TO 4                        AAUJ0547
      F3=FX                                                             AAUJ0548
      X3=X                                                              AAUJ0549
      U1=(F1-F2)/(X1-X2)                                                AAUJ0550
      U2=(F2-FX)/(X2-X)                                                 AAUJ0551
      CA=U1-U2                                                          AAUJ0552
      CB=(X1+X2)*U2-(X2+X)*U1                                           AAUJ0553
      CC=(X1-X)*F1-X1*(CA*X1+CB)                                        AAUJ0554
      IF(CA .EQ. 0.0) GO TO 8                                           AAUJ0555
      U3=0.5*CB/CA                                                      AAUJ0556
      U4=U3**2-CC/CA                                                    AAUJ0557
      IF(U4 .LT. 0.0) GO TO 4                                           AAUJ0558
      U5=SQRT(U4)                                                       AAUJ0559
      IF(X .GE. -U3) GO TO 10                                           AAUJ0560
      X=-U3-U5                                                          AAUJ0561
      GO TO 9                                                           AAUJ0562
   10 X=-U3+U5                                                          AAUJ0563
      GO TO 9                                                           AAUJ0564
    8 X=-CC/CB                                                          AAUJ0565
    9 IF(X .LT. XA) GO TO 4                                             AAUJ0566
      IF(X .GT. XB) GO TO 4                                             AAUJ0567
                                                                        AAUJ0568
C     TEST FOR OUTPUT                                                   AAUJ0569
                                                                        AAUJ0570
      R=ABS(X-X3)                                                       AAUJ0571
      R1=ABS(X-X2)                                                      AAUJ0572
      IF(R .GT. R1) R=R1                                                AAUJ0573
      EE=ABS(X)+E                                                       AAUJ0574
      IF(R/EE .GT. EPSI) GO TO 5                                        AAUJ0575
      MC=MC+1                                                           AAUJ0576
      G=FCN(X,RKA,BE2)                                                  AAUJ0577
      IF(MC .GT. MAXFUN) GO TO 17                                       AAUJ0578
      FX=G                                                              AAUJ0579
      IF(FX .EQ. 0.0) GO TO 18                                          AAUJ0580
      IF(FX*FA .LT. 0.0) GO TO 7                                        AAUJ0581
      XX=X+EPSI*EE                                                      AAUJ0582
      IF(XX .GE. XB) GO TO 18                                           AAUJ0583
      MC=MC+1                                                           AAUJ0584
      G=FCN(X,RKA,BE2)                                                  AAUJ0585
      IF(MC .GT. MAXFUN) GO TO 17                                       AAUJ0586
      FF=FCN(XX,RKA,BE2)                                                AAUJ0587
      FA=FF                                                             AAUJ0588
      XA=XX                                                             AAUJ0589
      GO TO 6                                                           AAUJ0590
    7 XX=X-EPSI*EE                                                      AAUJ0591
      IF(XX .LE. XA) GO TO 18                                           AAUJ0592
      MC=MC+1                                                           AAUJ0593
      FX=G                                                              AAUJ0594
      IF(MC .GT. MAXFUN) GO TO 17                                       AAUJ0595
      FF=FCN(XX,RKA,BE2)                                                AAUJ0596
      FB=FF                                                             AAUJ0597
      XB=XX                                                             AAUJ0598
    6 IF(FX*FF .GT. 0.0) GO TO 14                                       AAUJ0599
   18 R=EPSI*EE                                                         AAUJ0600
      FF=FCN(X,RKA,BE2)                                                 AAUJ0601
      RETURN                                                            AAUJ0602
   14 F1=F2                                                             AAUJ0603
      X1=X2                                                             AAUJ0604
      F2=FX                                                             AAUJ0605
      X2=X                                                              AAUJ0606
      X=XX                                                              AAUJ0607
      FX=FF                                                             AAUJ0608
      GO TO 3                                                           AAUJ0609
                                                                        AAUJ0610
    5 F1=F2                                                             AAUJ0611
      X1=X2                                                             AAUJ0612
      F2=F3                                                             AAUJ0613
      X2=X3                                                             AAUJ0614
      GO TO 1                                                           AAUJ0615
                                                                        AAUJ0616
   16 Continue ! WRITE(LU,301)
      RETURN                                                            AAUJ0618
                                                                        AAUJ0619
   17 Continue ! WRITE(LU,300) X,G,LFCN
  301 FORMAT(/10X,41H RZERO   FCN(A) AND FCN(B) HAVE SAME SIGN,/)       AAUJ0621
  300 FORMAT(/10X,38H RZERO   NUMBER OF ITERATIONS EXCEEDED,/           AAUJ0622
     110X,3H X=,E15.5,2X,8H FCN(X)=,E15.5,2X,6H LFCN=,I2/)              AAUJ0623
      RETURN                                                            AAUJ0624
                                                                        AAUJ0625
      END                                                               AAUJ0626

c-----------------------------------------------------------------------
