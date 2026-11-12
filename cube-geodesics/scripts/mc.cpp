// Brute-force Monte Carlo: geodesic statistics on the unit cube via min over all candidate unfoldings.
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <thread>
#include <vector>
struct Rng { uint64_t s[4];
  static uint64_t sm(uint64_t &x){ uint64_t z=(x+=0x9e3779b97f4a7c15ULL); z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL; z=(z^(z>>27))*0x94d049bb133111ebULL; return z^(z>>31);}
  Rng(uint64_t seed){ for(int i=0;i<4;i++) s[i]=sm(seed);}
  inline uint64_t next(){ uint64_t r=((s[1]*5)<<7|(s[1]*5)>>57)*9, t=s[1]<<17; s[2]^=s[0]; s[3]^=s[1]; s[1]^=s[2]; s[0]^=s[3]; s[2]^=t; s[3]=(s[3]<<45)|(s[3]>>19); return r;}
  inline double u(){ return ((next()>>11)+0.5)*(1.0/9007199254740992.0);} };
struct Acc { double n=0, n4=0, sopp=0, sopp2=0, nadj3=0, sadj=0, sadj2=0, ssame=0, ssame2=0; };
void work(uint64_t seed, uint64_t N, Acc* out){
  Rng r(seed); Acc A;
  for(uint64_t i=0;i<N;i++){
    double a=r.u(), b=r.u(), c=r.u(), d=r.u();
    double sq=[](double x){return x*x;}(0);
    #define S(x) ((x)*(x))
    double m3=std::min(std::min(S(b-d)+S(3-a-c), S(b-d)+S(1+a+c)), std::min(S(a-c)+S(3-b-d), S(a-c)+S(1+b+d)));
    double m4=std::min(std::min(std::min(S(3-d-a)+S(2-c-b), S(3-c-b)+S(2-d-a)), std::min(S(2-d+a)+S(1+c-b), S(2+c-b)+S(1-d+a))),
                       std::min(std::min(S(2+d-a)+S(1-c+b), S(2-c+b)+S(1+d-a)), std::min(S(1+d+a)+S(c+b), S(1+c+b)+S(d+a))));
    double dopp=std::sqrt(std::min(m3,m4)); A.n4+= (m4<m3); A.sopp+=dopp; A.sopp2+=dopp*dopp;
    // adjacent: reuse (a,b) as S on bottom, (c,d)=(y,z) on face x=1
    double d0=S(1+d-a)+S(c-b), d1=S(2-c-a)+S(1+d-b), d2=S(1+c-a)+S(d+b);
    double m=std::min(d1,d2); A.nadj3+=(m<d0); double da=std::sqrt(std::min(m,d0)); A.sadj+=da; A.sadj2+=da*da;
    double ds=std::sqrt(S(a-c)+S(b-d)); A.ssame+=ds; A.ssame2+=ds*ds;
    (void)sq;
  }
  A.n=N; *out=A;
}
int main(int argc,char**argv){
  uint64_t N=argc>1? strtoull(argv[1],0,10):1000000000ULL; int T=4;
  std::vector<Acc> acc(T); std::vector<std::thread> th;
  for(int t=0;t<T;t++) th.emplace_back(work, 12345+1000003ULL*t, N/T, &acc[t]);
  for(auto&x:th) x.join();
  Acc A; for(auto&x:acc){A.n+=x.n;A.n4+=x.n4;A.sopp+=x.sopp;A.sopp2+=x.sopp2;A.nadj3+=x.nadj3;A.sadj+=x.sadj;A.sadj2+=x.sadj2;A.ssame+=x.ssame;A.ssame2+=x.ssame2;}
  double p4=A.n4/A.n, pa=A.nadj3/A.n;
  auto se=[&](double s,double s2){double m=s/A.n; return std::sqrt((s2/A.n-m*m)/A.n);};
  printf("N=%.0f\n",A.n);
  printf("P4 (opposite, four-face)   = %.8f +- %.8f\n",p4,std::sqrt(p4*(1-p4)/A.n));
  printf("P  (adjacent, third face)  = %.8f +- %.8f\n",pa,std::sqrt(pa*(1-pa)/A.n));
  printf("mean geodesic opposite     = %.8f +- %.8f\n",A.sopp/A.n,se(A.sopp,A.sopp2));
  printf("mean geodesic adjacent     = %.8f +- %.8f\n",A.sadj/A.n,se(A.sadj,A.sadj2));
  printf("mean same face             = %.8f +- %.8f\n",A.ssame/A.n,se(A.ssame,A.ssame2));
  double D=(A.ssame/A.n)/6+4*(A.sadj/A.n)/6+(A.sopp/A.n)/6;
  printf("mean geodesic whole surface= %.8f\n",D);
}
