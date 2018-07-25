/*
This main asks for the Poisson parameter lambda, then generates 10^8
random Poisson-lambda variates using the compacted 5-table method.
It then does a goodness-of-fit test on the output.
It then asks for the binomial parameters n,p, generates 10^8
binomial(n,p) variates with the compacted 5-table method and tests
the output.
Finally, it asks for the hypergeometric parameters, then generates and
does a test on the output of 10^8 such hypergeometric variates.
Generation speed is 10-15 nanoseconds for each variate,
some 70-100 million per second.
Output is displayed on the screen and also sent to the file tests.out.
*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>

# define dg(m,k) ((m>>(30-6*k))&63)  /* gets kth digit of m (base 64) */

static int *P;   /* Probabilities as an array of 30-bit integers*/
static int size,t1,t2,t3,t4,offset,last;  /* size of P[], limits for table lookups */
static unsigned long jxr=182736531; /* Xorshift RNG */
static short int *AA,*BB,*CC,*DD,*EE;      /* Tables for condensed table-lookup */
static FILE *fp;

void get5tbls() /* Creates the 5 tables after array P[size] has been created */
{ int i,j,k,m,na=0,nb=0,nc=0,nd=0,ne=0;
    /* get table sizes, malloc */
for(i=0;i<size;i++){m=P[i];na+=dg(m,1);nb+=dg(m,2);nc+=dg(m,3);nd+=dg(m,4);ne+=dg(m,5);}
AA=malloc(na*sizeof(int));BB=malloc(nb*sizeof(int));
CC=malloc(nc*sizeof(int));DD=malloc(nd*sizeof(int));EE=malloc(ne*sizeof(int));
printf(" Table sizes:%d,%d,%d,%d,%d,total=%d\n",na,nb,nc,nd,ne,na+nb+nc+nd+ne);
fprintf(fp," Table sizes:%d,%d,%d,%d,%d,total=%d\n",na,nb,nc,nd,ne,na+nb+nc+nd+ne);
t1=na<<24; t2=t1+(nb<<18); t3=t2+(nc<<12); t4=t3+(nd<<6);
na=nb=nc=nd=ne=0;
   /* Fill tables AA,BB,CC,DD,EE */
for(i=0;i<size;i++){m=P[i]; k=i+offset;
         for(j=0;j<dg(m,1);j++) AA[na+j]=k; na+=dg(m,1);
         for(j=0;j<dg(m,2);j++) BB[nb+j]=k; nb+=dg(m,2);
         for(j=0;j<dg(m,3);j++) CC[nc+j]=k; nc+=dg(m,3);
         for(j=0;j<dg(m,4);j++) DD[nd+j]=k; nd+=dg(m,4);
         for(j=0;j<dg(m,5);j++) EE[ne+j]=k; ne+=dg(m,5);
                   }
}/* end get5tbls */

void PoissP(double lam)  /* Creates Poisson Probabilites */
{int i,j=-1,nlam;        /* P's are 30-bit integers, assumed denominator 2^30 */
double p=1,c,t=1.;
    /* generate  P's from 0 if lam<21.4 */
if(lam<21.4){p=t=exp(-lam); for(i=1;t*2147483648.>1;i++) t*=(lam/i);
          size=i-1; last=i-2;
          /* Given size, malloc and fill P array, (30-bit integers) */
          P=malloc(size*sizeof(int)); P[0]=exp(-lam)*(1<<30)+.5;
          for(i=1;i<size;i++) {p*=(lam/i); P[i]=p*(1<<30)+.5; }
            }
    /* If lam>=21.4, generate from largest P up,then largest down */
if(lam>=21.4)
{nlam=lam;      /*first find size */
c=lam*exp(-lam/nlam);
for(i=1;i<=nlam;i++) t*=(c/i);
p=t;
for(i=nlam+1;t*2147483648.>1;i++) t*=(lam/i);
last=i-2;
t=p; j=-1;
for(i=nlam-1;i>=0;i--){t*=((i+1)/lam); if(t*2147483648.<1){j=i;break;} }
offset=j+1;  size=last-offset+1;
   /* malloc and fill P array as 30-bit integers */
P=malloc(size*sizeof(int));
t=p; P[nlam-offset]=p*(1<<30)+0.5;
for(i=nlam+1;i<=last;i++){t*=(lam/i); P[i-offset]=t*(64<<24)+.5;}
t=p;
for(i=nlam-1;i>=offset;i--){t*=((i+1)/lam);P[i-offset]=t*(1<<30)+0.5;}
} /* end lam>=21.4*/
} /*end PoissP*/

void BinomP(int n, double p)  /*Creates Binomial Probabilities */
         /* Note: if n large and p near 1, generate j=Binom(n,1-p), return n-j*/
{double h,t,pm,q;
 int i,j,m;
/* first find size of P array */
m=(n+1)*p;q=1-p;h=p/q;
for(i=0,j=0,pm=1;!(i==n-m&&j==m);)
{if (pm<1&&i<n-m) {pm *= q*(n-i)/(n-m-i);i++;}
else {pm *= p;j++;}}
for(i=m+1,t=pm;t*2147483648.>1;i++)
t *= h*(n-i+1)/i;
last=i-2;
for(i=m-1,t=pm,j=-1;i>=0;i--)
{t *= (i+1)/h/(n-i);
if(t*2147483648.<1) {j=i;break;} }
offset=j+1;
size=last-offset+1;
/* Then malloc and assign P values as 30-bit integers */
P=(int *)malloc(size*sizeof(int));
P[m-offset]=pm*(1<<30)+0.5;
for(i=m+1,t=pm;i<=last;i++){t *= h*(n-i+1)/i; P[i-offset]=t*(64<<24)+.5;}
for(i=m-1,t=pm;i>=offset;i--){t /= h*(n-i)/(i+1);P[i-offset]=t*(1<<30)+0.5;}
} /*end BinomP*/

void HyperGeometricP (long int g, long int b, long int s){
long int i,j,k,mode;
double pmode,t,modef;
modef = (double)(g * s - b + s - 1) / (g + b + 2.);
mode = (long int)modef + 1;

for (i=0,j=0,k=0,pmode=1; !(i==mode && j==s-mode && k==s);){
if (pmode > 1 && k<s){pmode /= (double)(b+g-s+k+1)/(k+1);k++;}
else if (i<mode){pmode *= (double)(g-mode+i+1)/(i+1);i++;}
	 else if (j<s-mode){pmode *= (double)(b-s+mode+j+1)/(j+1);j++;}
	 else if (k<s){pmode /= (double)(b+g-s+k+1)/(k+1);k++;} }

t=pmode;
for(i=mode+1;t*2147483648.>1;i++) t *= ((double)(g-i+1)*(s-i+1)/i/(b-s+i));
last=i-2;
t=pmode;
j=-1;
for(i=mode-1;i>=0;i--)
{t /= ((double)(g-i)*(s-i)/(i+1)/(b-s+i+1));
 if(t*2147483648.<1) {j=i;break;} }
offset=j+1;
size=last-offset+1;

P=(int *)malloc(size*sizeof(int));
t=pmode; P[mode-offset]=pmode*(1<<30)+0.5;
for(i=mode+1;i<=last;i++){t *= ((double)(g-i+1)*(s-i+1)/i/(b-s+i));
P[i-offset]=t*(64<<24)+.5;}
t=pmode;
for(i=mode-1;i>=offset;i--){t /= ((double)(g-i)*(s-i)/(i+1)/(b-s+i+1));
P[i-offset]=t*(1<<30)+0.5;}
}





/* Discrete random variable generating function */

int Dran() /* Uses 5 compact tables */
{unsigned long j;
 jxr^=jxr<<13; jxr^=jxr>>17; jxr^=jxr<<5; j=(jxr>>2);
if(j<t1) return AA[j>>24];
if(j<t2) return BB[(j-t1)>>18];
if(j<t3) return CC[(j-t2)>>12];
if(j<t4) return DD[(j-t3)>>6];
return EE[j-t4];
}


double Phi(double x)
 {long double s=x,t=0,b=x,x2=x*x,i=1;
    while(s!=t) s=(t=s)+(b*=x2/(i+=2));
    return  .5+s*exp(-.5*x2-.91893853320467274178L);
 }

double chisq(double z,int n)
{double s=0.,t=1.,q,h;
 int i;
 if(z<=0.) return (0.);
 if(n>3000) return Phi((exp(log(z/n)/3.)-1.+2./(9*n))/sqrt(2./(9*n)));
 h=.5*z;
 if(n&1){q=sqrt(z); t=2*exp(-.5*z-.918938533204673)/q;
    for(i=1;i<=(n-2);i+=2){ t=t*z/i; s+=t;}
    return(2*Phi(q)-1-s);
        }
 for(i=1;i<n/2;i++) { t=t*h/i; s+=t;}
 return (1.-(1+s)*exp(-h));
}


void Dtest(n)
/* requires static 'size', static int array P[size] */
/* generates n Dran()'s, tests output */
{ double x=0,y,s=0,*E;
  int kb=0,ke=1000,i,j=0,*M;
 E=malloc(size*sizeof(double));
 M=malloc(size*sizeof(int));
 for(i=0;i<size;i++) {E[i]=(n+0.)*P[i]/1073741824.;M[i]=0;}
 s=0; for(i=0;i<size;i++) {s+=E[i]; if(s>10){kb=i;E[kb]=s;break;} }
 s=0; for(i=size-1;i>0;i--) {s+=E[i]; if(s>10){ke=i;E[ke]=s;break;} }
 s=0; x=0; for(i=0;i<n;i++)
 {j=Dran(); if(j<kb+offset) j=kb+offset; if(j>ke+offset) j=ke+offset;
  M[j-offset]++;}

 printf("\n   D     Observed     Expected    (O-E)^2/E   sum\n");
 fprintf(fp,"\n   D     Observed     Expected    (O-E)^2/E   sum\n");
 for(i=kb;i<=ke;i++){ y=M[i]-E[i]; y=y*y/E[i]; s+=y;

 printf("%4d %10d  %12.2f   %7.2f   %7.2f\n",i+offset,M[i],E[i],y,s);
 fprintf(fp,"%4d %10d   %12.2f   %7.2f   %7.2f\n",i+offset,M[i],E[i],y,s);
                  }
 printf("    chisquare for %d d.f.=%7.2f, p=%7.5f\n",ke-kb,s,chisq(s,ke-kb));
 fprintf(fp,"     chisquare for %d d.f.=%7.2f, p=%7.5f\n",ke-kb,s,chisq(s,ke-kb));
}

int main(){
int j=0,n,g=0,b=0,s=1,nsmpls=100000000;
double lam,p;
fp=fopen("tests.out","w");
printf("  Enter lambda:\n");
scanf("%lf",&lam);
PoissP(lam); get5tbls();
Dtest(nsmpls);
fprintf(fp," Above results for sample of %d from Poisson, lambda=%3.2f\n\n",nsmpls,lam);
free(P);
printf("\n Enter n and p for Binomial:\n");
scanf("%d %lf",&n,&p);
BinomP(n,p);get5tbls();
Dtest(nsmpls);
fprintf(fp," Above result for sample of %d from Binomial(%d,%3.3f)\n\n",nsmpls,n,p);
free(P);
while(s>g||s>b){
printf("\n Enter N1, N2 and K (K<=N1,N2) for HyperGeometric:\n");
scanf("%d %d %d",&g,&b,&s);}
HyperGeometricP (g,b,s);get5tbls();
Dtest(nsmpls);
fprintf(fp," Above result for sample of %d from HyperGeometric(%d,%3.3f)\n\n",nsmpls,n,p);
free(P);
printf(" Test results sent to file tests.out\n");
return 0;
}


