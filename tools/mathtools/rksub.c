
/* RKSUB.C     Gerhard Heinzel   24.4.92






  Subroutinen zur numerischen Integration gewoehnlicher

  Differentialgleichungen mit Runge-Kutta-Tripeln, mit Mîglichkeit
   der Ausgabe vieler ("dichter") Zwischenergebnisse.
   Sprache : Microsoft C 6.0
   Beschreibung und Formeln in c't 8/92                            */

#include "rksub.h"

long ngood,nfail;    /* SchrittzÑhler: extern sichtbar in rktest.c */

/*                  Koeffizienten des Tripels, gesetzt von rkcoeff */
static int s,ss,q,p,d,fsal;        /* int-Konstanten aus Tabelle 3 */
static DOUBLE a[SMAX][SMAX-1],b[SMAX],bdach[SMAX],berr[SMAX],
    c[SMAX],beta[DMAX+1][SMAX], potenz;     /* Tripel-Koeffizenten */

/*                        Parameter, die mit rkinit gesetzt werden */
static DOUBLE eps,tend;           /* Genauigkeit, Integrationsende */
static int nvar,autonom;       /* Systemgrî·e, t-abhÑngig ja/nein? */
static void (*f)(DOUBLE, DOUBLE*, DOUBLE*);      /* Zeiger auf DGL */

/*                          FÅr Interpolation gesicherte Variablen */
static int neu;                       /* Polynom schon berechnet ? */
static DOUBLE k[SMAX][NVARMAX];              /* Zwischenergebnisse */
static DOUBLE h0;                           /* Letzte Schrittweite */
static DOUBLE tt;                              /* Letzte Startzeit */
static DOUBLE ybeg[NVARMAX];   /* Startpunkt des letzten Schrittes */

/*                          sonstige dauerhaft benîtigte Variablen */
static DOUBLE h;                              /* Neue Schrittweite */
static DOUBLE tlast=DBL_MAX;      /* Endzeit des letzten Schrittes */


void rkcoeff(void)          /* Koeffizienten des Tripels einsetzen */
{
register int i,j;

/*                           Koeffizienten auf Null initialisieren */
for (i=0; i<SMAX; i++) {
    b[i]=bdach[i]=berr[i]=c[i]=0.0;
    for (j=0; j<SMAX-1; j++)
        a[i][j]=0.0;
    for (j=0; j<DMAX+1; j++)
        beta[j][i]=0.0;
    }

/***************** Anfang auszutauschender Teil ********************/
/* DP8(6)T7 Runge-Kutta-Tripel von Dormand/Prince
   Koeffizienten zur VerfÅgung gestellt von Dr. John Dormand
   umgeformt (30 Stellen) von Gerhard Heinzel MÑrz 92
   alle nicht genannten Koeffizienten sind Null! */

q=8; p=6; s=12; ss=16; d=7; fsal=NO;

a[1][0]=0.0526001519587677318785587544488L;
a[2][0]=0.0197250569845378994544595329183L;
a[2][1]=0.0591751709536136983633785987549L;
a[3][0]=0.0295875854768068491816892993775L;
a[3][2]=0.0887627564304205475450678981324L;
a[4][0]=0.241365134159266685502369798665L;
a[4][2]=-0.884549479328286085344864962717L;
a[4][3]=0.924834003261792003115737966543L;
a[5][0]=0.037037037037037037037037037037L;
a[5][3]=0.170828608729473871279604482173L;
a[5][4]=0.125467687566822425016691814123L;
a[6][0]=0.037109375L;
a[6][3]=0.170252211019544039314978060272L;
a[6][4]=0.0602165389804559606850219397283L;
a[6][5]=-0.017578125L;
a[7][0]=0.0370920001185047927108779319836L;
a[7][3]=0.170383925712239993810214054705L;
a[7][4]=0.107262030446373284651809199168L;
a[7][5]=-0.0153194377486244017527936158236L;
a[7][6]=0.00827378916381402288758473766002L;
a[8][0]=0.624110958716075717114429577812L;
a[8][3]=-3.36089262944694129406857109825L;
a[8][4]=-0.868219346841726006818189891453L;
a[8][5]=27.5920996994467083049415600797L;
a[8][6]=20.1540675504778934086186788979L;
a[8][7]=-43.4898841810699588477366255144L;
a[9][0]=0.477662536438264365890433908527L;
a[9][3]=-2.48811461997166764192642586468L;
a[9][4]=-0.590290826836842996371446475743L;
a[9][5]=21.2300514481811942347288949897L;
a[9][6]=15.2792336328824235832596922938L;
a[9][7]=-33.2882109689848629194453265587L;
a[9][8]=-0.0203312017085086261358222928593L;
a[10][0]=-0.93714243008598732571704021658L;
a[10][3]=5.18637242884406370830023853209L;
a[10][4]=1.09143734899672957818500254654L;
a[10][5]=-8.14978701074692612513997267357L;
a[10][6]=-18.5200656599969598641566180701L;
a[10][7]=22.7394870993505042818970056734L;
a[10][8]=2.49360555267965238987089396762L;
a[10][9]=-3.0467644718982195003823669022L;
a[11][0]=2.27331014751653820792359768449L;
a[11][3]=-10.5344954667372501984066689879L;
a[11][4]=-2.00087205822486249909675718444L;
a[11][5]=-17.9589318631187989172765950534L;
a[11][6]=27.9488845294199600508499808837L;
a[11][7]=-2.85899827713502369474065508674L;
a[11][8]=-8.87285693353062954433549289258L;
a[11][9]=12.3605671757943030647266201528L;
a[11][10]=0.643392746015763530355970484046L;
a[12][0]=0.0542937341165687622380535766363L;
a[12][5]=4.45031289275240888144113950566L;
a[12][6]=1.89151789931450038304281599044L;
a[12][7]=-5.8012039600105847814672114227L;
a[12][8]=0.31116436695781989440891606237L;
a[12][9]=-0.152160949662516078556178806805L;
a[12][10]=0.201365400804030348374776537501L;
a[12][11]=0.0447106157277725905176885569043L;
a[13][0]=0.0561675022830479523392909219681L;
a[13][6]=0.253500210216624811088794765333L;
a[13][7]=-0.246239037470802489917441475441L;
a[13][8]=-0.124191423263816360469010140626L;
a[13][9]=0.15329179827876569731206322685L;
a[13][10]=0.00820105229563468988491666602057L;
a[13][11]=0.00756789766054569976138603589584L;
a[13][12]=-0.008298L;
a[14][0]=0.0318346481635021405060768473261L;
a[14][5]=0.0283009096723667755288322961402L;
a[14][6]=0.0535419883074385676223797384372L;
a[14][7]=-0.0549237485713909884646569340306L;
a[14][10]=-0.000108347328697249322858509316994L;
a[14][11]=0.000382571090835658412954920192323L;
a[14][12]=-0.000340465008687404560802977114492L;
a[14][13]=0.141312443674632500278074618366L;
a[15][0]=-0.428896301583791923408573538692L;
a[15][5]=-4.69762141536116384314449447206L;
a[15][6]=7.68342119606259904184240953878L;
a[15][7]=4.06898981839711007970213554331L;
a[15][8]=0.356727187455281109270669543021L;
a[15][12]=-0.00139902416515901462129418009734L;
a[15][13]=2.9475147891527723389556272149L;
a[15][14]=-9.15095847217987001081870187138L;
b[0]=0.0542937341165687622380535766363L;
b[5]=4.45031289275240888144113950566L;
b[6]=1.89151789931450038304281599044L;
b[7]=-5.8012039600105847814672114227L;
b[8]=0.31116436695781989440891606237L;
b[9]=-0.152160949662516078556178806805L;
b[10]=0.201365400804030348374776537501L;
b[11]=0.0447106157277725905176885569043L;
bdach[0]=0.0633186814225382321037030864969L;
bdach[5]=2.L;
bdach[6]=1.13526283625113359349971992899L;
bdach[7]=-2.71535512663407182042077739032L;
bdach[8]=0.0540072979937810585557526143181L;
bdach[9]=0.206914334754189019481881265971L;
bdach[10]=0.211141360484657326262031937633L;
bdach[11]=0.0447106157277725905176885569043L;
beta[0][0]=1.L;
beta[1][0]=-10.2660570737593065784211883858L;
beta[2][0]=48.1618509685664566301953957035L;
beta[3][0]=-114.933048749978332538237198113L;
beta[4][0]=147.464468756697683076313870249L;
beta[5][0]=-97.0668536301136808309254120056L;
beta[6][0]=25.6939334627037490033125861294L;
beta[1][5]=13.9176536317766044139487363529L;
beta[2][5]=-154.787872666637155968890178909L;
beta[3][5]=522.921908960821874913658927437L;
beta[4][5]=-456.259188402087812547255490252L;
beta[5][5]=-75.5319373213575356705607913937L;
beta[6][5]=154.189748690236433740539936271L;
beta[1][6]=2.60560375199360945784871749873L;
beta[2][6]=-21.6228223846265042267806054463L;
beta[3][6]=2.5351820289667551481771305393L;
beta[4][6]=292.254174659904062526209972298L;
beta[5][6]=-505.409999332968918197772789988L;
beta[6][6]=231.529379176045495675360391089L;
beta[1][7]=-15.0189442235196845156288192986L;
beta[2][7]=160.094477089730476115815838246L;
beta[3][7]=-474.307182603764347814837943467L;
beta[4][7]=135.960369161738372873086881756L;
beta[5][7]=545.109194526418722342950330441L;
beta[6][7]=-357.6391179106141237828534991L;
beta[1][8]=3.05052768331848795994226318181L;
beta[2][8]=-38.5439672918906325046616718257L;
beta[3][8]=174.471400092198840731590697256L;
beta[4][8]=-337.051347023877126426419627091L;
beta[5][8]=291.789875090832560137864946245L;
beta[6][8]=-93.4053241836243100039076917036L;
beta[1][9]=-1.32787443276552122773643549547L;
beta[2][9]=16.6617704300495419971752379106L;
beta[3][9]=-74.4402781412630338777638949861L;
beta[4][9]=140.752100161916063360485884733L;
beta[5][9]=-119.25620210405119948759211032L;
beta[6][9]=37.4583231364516331568751393513L;
beta[1][10]=2.84453363267287932097850677632L;
beta[2][10]=-36.5582954899101192708375071275L;
beta[3][10]=170.690071691475136612404224555L;
beta[4][10]=-345.974848548049551057433757452L;
beta[5][10]=313.29955362357798519473577163L;
beta[6][10]=-104.099649508962300451472461844L;
beta[1][11]=0.765710625952786589708768163779L;
beta[2][11]=-9.90699553561936636937481161302L;
beta[3][11]=46.8029919188743947246274324406L;
beta[4][11]=-96.5198694669957042802037349345L;
beta[5][11]=88.7431665001761650491043980787L;
beta[6][11]=-29.8402934266605031233443635787L;
beta[1][12]=-1.08899033645133331082069811682L;
beta[2][12]=14.0970130423200021011792868748L;
beta[3][12]=-66.6823059129436396177354540236L;
beta[4][12]=137.962990634743749929648014948L;
beta[5][12]=-127.822164017679922856703324741L;
beta[6][12]=43.5334565900111437544321750583L;
beta[1][13]=18.1485055208547272566564049647L;
beta[2][13]=-127.633109492538752948863041185L;
beta[3][13]=357.341951612965727834419198923L;
beta[4][13]=-500.70315079092238879726984475L;
beta[5][13]=349.170357108828969603452232647L;
beta[6][13]=-96.3245539591882829483949506004L;
beta[1][14]=-9.19463239247835540004519844455L;
beta[2][14]=93.3567459327893934316789144061L;
beta[3][14]=-282.627261870436320846613623435L;
beta[4][14]=361.140077188033322163602783601L;
beta[5][14]=-201.852190533523478513854362299L;
beta[6][14]=39.1772616756154391652314861716L;
beta[1][15]=-4.43603638759489396643105719702L;
beta[2][15]=56.6812053977666610133631429662L;
beta[3][15]=-261.773429026917055269689497125L;
beta[4][15]=520.974223668899329179235046895L;
beta[5][15]=-461.172799910139666770698888295L;
beta[6][15]=149.726836257985625814221252755L;

/********************** Ende spezieller Teil **********************/

for (i=0; i<ss; i++)
    for (j=0; j<i; j++)
        c[i] += a[i][j];                           /* Formel (15) */
for (i=0; i<s; i++)
    berr[i]=bdach[i]-b[i];
potenz = 1.L/(p+1.L);                  /* Exponent in Formel (21) */
}


void rkinit(void (*f_p)(DOUBLE, DOUBLE*, DOUBLE*),int autonom_p,
            int nvar_p, DOUBLE eps_p, DOUBLE tend_p, DOUBLE h_p)
{
/*                        öbergebene Parameter als static sichern */
eps=eps_p; nvar=nvar_p; tend=tend_p; autonom=autonom_p;
f=f_p;                                     /* Zeiger auf Funktion */


if (FABS(h_p) >= 2*EPSILON)               /* falls ungleich Null, */
    h=h_p;                       /* Startwert Schrittweite setzen */
}


DOUBLE rkstep(DOUBLE t, DOUBLE y[NVARMAX], DOUBLE yneu[NVARMAX])
{
DOUBLE ytemp[NVARMAX],ydiff,err,hfakt,ttemp=0.L;
register int i,j,l;
int fail,ifail=0;

tt=t;                     /* Anfangszeit fÅr Interpolation sichern */
if (!fsal || FABS(t-tlast)>2.*EPSILON)
    (*f)(t,y,k[0]);                            /* Erste Auswertung */
else                                  /* FSAL : k0 schon berechnet */
    memcpy(k[0],k[s-1],nvar*sizeof(DOUBLE));
do {                         /* Schritt ausfÅhren, bis genau genug */
    if (t+h>tend)                             /* Letzter Schritt ? */
        h=tend-t;                            /* Dann genau treffen */
    for (i=1; i<s; i++) {  /* k1...ks-1 berechnen, besser explizit */
        for (l=0; l<nvar; l++) {      /* l zÑhlt Vektorkomponenten */
          ytemp[l]=0.L;
          for (j=0; j<i; j++)                  /* y-Argument fÅr f */
            ytemp[l] += a[i][j]*k[j][l];
          ytemp[l] = y[l] + h * ytemp[l];           /* Formel (14) */
          }
        if (!autonom)                                /* t-Argument */
          ttemp=t+c[i]*h;
        (*f)(ttemp,ytemp,k[i]);            /* Funktion f auswerten */
        }
    err=0.L;    /* Erbebnis und FehlerabschÑtzung dieses Schrittes */
    for (l=0; l<nvar; l++) {          /* l zÑhlt Vektorkomponenten */
        yneu[l] = ydiff = 0.L;
        for (i=0; i<s; i++) {
          yneu[l] += b[i]*k[i][l];       /* Zwischensumme Ergebnis */
          ydiff += berr[i]*k[i][l];      /* dito FehlerabschÑtzung */
          }
        ydiff = FABS(h*ydiff);
        err = max(err,ydiff);                     /* Maximums-Norm */
        yneu[l] = y[l] + h * yneu[l];                  /* Ergebnis */
        }
    if ((fail=(err>eps)) == YES)      /* Schritt zu grob gewesen ? */
        ++ifail;                     /* ZÑhler erfolglose Schritte */
    if (ifail>5)
        printf("!");

    hfakt = 0.9L * POW(eps/err,potenz);             /* Formel (21) */
    if (hfakt < 0.1L)           /* SchrittweitenÑnderung begrenzen */
        hfakt = 0.1L;
    if (hfakt > 5.L)
        hfakt = 5.L;
    h0=h;           /* Alte Schrittweite fÅr Interpolation sichern */
    if (!ifail || hfakt < 1.0L)           /* Schrittweite anpassen */
        h *= hfakt;
} while(fail);
nfail+=ifail; ++ngood;                                /* Statistik */
memcpy(ybeg,y,nvar*sizeof(DOUBLE)); /* y sichern fÅr Interpolation */
neu=YES;         /* Interpolationspolynom mu· neu berechnet werden */
return (tlast = t+h0);   /* Zeit am Ende des Schrittes zurÅckgeben */
}

void rkipol(DOUBLE ts, DOUBLE ys[NVARMAX])
{
DOUBLE sigma,ytemp[NVARMAX],ttemp;
static DOUBLE aa[DMAX+1][NVARMAX];   /* Koeffizienten des Polynoms */
register int i,j,l;

if (neu) {                     /* Polynom mu· neu berechnet werden */
    for (i=s; i<ss; i++) { /* ks...ks-1 berechnen, besser explizit */
        for (l=0; l<nvar; l++) {    /* Sparmîglichkeit fÅr i==ss-1 */
          ytemp[l]=0.L;
          for (j=0; j<i; j++)
            ytemp[l] += a[i][j]*k[j][l];
          ytemp[l] = ybeg[l] + h0 * ytemp[l];       /* Formel (24) */
          }
        if (!autonom)
          ttemp=tt+c[i]*h0;
        (*f)(ttemp,ytemp,k[i]);              /* Funktion auswerten */
        }
/*                                 Polynom-Koeffizienten berechnen */
    for (l=0; l<nvar; l++) {         /* Sparmîglichkeit siehe Text */
        aa[0][l] = ybeg[l];                         /* Formel (29) */
        for (j=1; j<=d; j++) {
          aa[j][l]=0.L;
          for (i=0; i<ss; i++)
            aa[j][l] += beta[j-1][i]*k[i][l];       /* Formel (30) */
          aa[j][l] *= h0;
          }
        }
    neu=NO;
    }
sigma = (ts-tt)/h0;                                 /* Formel (25) */
if (sigma<0.L || sigma>1.L)
    printf("sigma Bereich 0...1");
for (l=0; l<nvar; l++) {             /* Sparmîglichkeit siehe Text */
    ys[l]=aa[d][l];                 /* Formel (28) : Horner-Schema */
    for (i=d-1; i>=0; i--)
      ys[l]=ys[l]*sigma+aa[i][l];
    }
}























