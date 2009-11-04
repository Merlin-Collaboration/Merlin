#ifndef _collimatortable
#define _collimatortable

#include <iostream>
#include <fstream>
using namespace std;

class collimatortable{
float* coeff;
float step; // step size in s
int ncoeff; // number 
float lo,hi;
public:

~collimatortable(){ delete [] coeff;};

bool inrange(double x){return (x>=lo) && (x<= hi);}

double parabolic(double v1, double v2, double v3,double d){
     return v2+d*(v3-v1)/2+d*d*(v3+v1-2*v2)/2;
   }

double interpolate(double s){
 int index=int(0.5+s/step);
 if(index<0) index=0;
 int i1;
 if(index==0){i1=1;} else if(index==ncoeff-1) {i1=ncoeff-2;} else {i1=index;}
 return(index>=ncoeff)?0:parabolic(coeff[i1-1],coeff[i1],coeff[i1+1],s/step-i1);
 }

collimatortable(char* file, double Gamma=0, double xi=0){
float lo1,hi1,lo3,hi3;
ifstream f;
f.open(file);
if(!f) {cout<<" cannot open file"<<file<<endl; exit(1);}
int n1,n2,n3;
f>>n3>>lo3>>hi3;
float step3=(hi3-lo3)/(n3-1);
f>>n2>>lo>>hi;
step=(hi-lo)/(n2-1);
f>>n1>>lo1>>hi1;
float step1=(hi1-lo1)/(n1-1);
ncoeff=n2;
coeff=new float[n2];
float array[n1][n2][n3];

while(f.get()!='{') ;
for(int i1=0;i1<n1;i1++){
   while(f.get()!='{') ;
   for(int i2=0;i2<n2;i2++){
     while(f.get()!='{') ;
     for(int i3=0;i3<n3;i3++){
//         cout<<i1<<" "<<i2<<" "<<i3<<flush;
         f>>array[i1][i2][i3];
  //       cout<<" "<<f<<endl;
         if(i3<n3-1) while(f.get()!=',');
         }
     while(f.get()!='}') ;
     if(i2<n2-1) while(f.get()!=',');
     }
   while(f.get()!='}') ;
   if(i1<n1-1) while(f.get()!=',');
   }
   while(f.get()!='}') ;
  // cout<<" Table read OK "<<endl;

if((Gamma==0) && (xi==0)) { // In principle the inerpolation should 
                            // handle this but it seems cleaner this way
   for(int i=0;i<ncoeff;i++) {coeff[i]=array[0][i][0];}
   }
   else if(Gamma==0) { // interpolate xi value
     int index=int(0.5+xi/step1);
     if(index<0) index=0;
     int i1;
     if(index==0){i1=1;} else if(index==n1-1) {i1=n1-2;} else {i1=index;}
     for(int i=0; i<ncoeff; i++){
     if(index>=n1) {coeff[i]=0;} else 
      {
      double d=(xi-i1*step1)/step1;
      coeff[i]=parabolic(array[0][i][i1-1],array[0][i][i1],array[0][i][i1+1],d);
     }
     }
}
else if(xi==0) { // interpolate Gamma value 
     int index=int(0.5+Gamma/step3);
     if(index<0) index=0;
     int i1;
     if(index==0){i1=1;} else if(index==n3-1) {i1=n3-2;} else {i1=index;}
     for(int i=0; i<ncoeff; i++){
     if(index>=n1) {coeff[i]=0;} else 
      {
      double d=(Gamma-i1*step3)/step3;
      coeff[i]=parabolic(array[i1-1][i][0],array[i1][i][0],array[i1+1][i][0],d);
     }
     }
}
else{// 2D interpolation
     int indg=int(0.5+Gamma/step3);
     if(indg<0) indg=0;
     int ig1;
     if(indg==0){ig1=1;} else if(indg==n1-1) {ig1=n1-2;} else {ig1=indg;}
     int indc=int(0.5+xi/step1);
     if(indc<0) indc=0;
     int ic1;
     if(indc==0){ic1=1;} else if(indc==n3-1) {ic1=n3-2;} else {ic1=indc;}
      double dg=(Gamma-ig1*step3)/step3;
      double dc=(xi-ic1*step1)/step1;
      for(int i=0; i<ncoeff; i++){
       float f0= array[ig1][i][ic1];
       float fc=(array[ig1][i][ic1+1]-array[ig1][i][ic1-1])/2;
       float fg=(array[ig1+1][i][ic1]-array[ig1-1][i][ic1])/2;
       float fcc=(array[ig1][i][ic1+1]+array[ig1][i][ic1-1])/2-array[ig1][i][ic1];
       float fgg=(array[ig1+1][i][ic1]+array[ig1-1][i][ic1])/2-array[ig1][i][ic1];
       float apm=array[ig1+1][i][ic1-1]-(f0+fg+fgg-fc+fcc);
       float app=array[ig1+1][i][ic1+1]-(f0+fg+fgg+fc+fcc);
       float amm=array[ig1-1][i][ic1-1]-(f0-fg+fgg-fc+fcc);
       float amp=array[ig1-1][i][ic1+1]-(f0-fgg+fc+fcc);
       float fgc=(app+amm-apm-amp)/4.;
       float fggc=(app-amp+apm-amm)/4.;
       float fgcc=(app+amp-apm-amm)/4.;
       float fggcc=(apm+app+amp+amm)/4.;
   //  cout<<" check "<<f0<<" "<<array[ig1][i][ic1]<<endl;
   //  cout<<" check "<<f0+fg+fgg<<" "<<array[ig1+1][i][ic1]<<endl;
  //   cout<<" check "<<f0-fg+fgg<<" "<<array[ig1-1][i][ic1]<<endl;
  //   cout<<" check "<<f0-fc+fcc<<" "<<array[ig1][i][ic1-1]<<endl;
  //   cout<<" check "<<f0+fc+fcc<<" "<<array[ig1][i][ic1+1]<<endl;
  //   cout<<" check "<<f0+fc+fcc+fg+fgg+fgc+fggc+fgcc+fggcc<<" "<<array[ig1+1][i][ic1+1]<<endl;
  //   cout<<" check "<<f0-fc+fcc-fg+fgg+fgc-fggc-fgcc+fggcc<<" "<<array[ig1-1][i][ic1-1]<<endl;
    
     coeff[i]=f0+dg*fg+dc*fc+dg*dg*fgg+dc*dc*fcc+dg*dc*fgc+dg*dg*dc*fggc
     +dg*dc*dc*fgcc+dg*dg*dc*dc*fggcc;
     }
}

cout<<" File "<<file<<" read and table calculated for Gamma of "<<Gamma<<" and xi of "<<xi<<endl;
}
};

#endif
