#include "mathe.h"

/* returns a float between 0 and 1 (inclusive)*/
float myrand(){
  return (float)(((double)rand())/RAND_MAX);
}
///* computes the inverse of a 3x3 matrix A, where
// *   [a[0] a[3] a[6]]
// * A=[a[1] a[4] a[7]]
// *   [a[2] a[5] a[8]]
// * stores the inverse in dest if dest!=NULL
// * stores the determinant in determinant if non-null
// * returns the inverseof A if A invertible, NULL otherwise
// * Warning: check the determinant to make sure the inverse
// * is OK (Cf determinant close to 0, but not exactly equal
// * to 0)*/
//float* invert3x3(float *a, float *dest, float *determinant){
//  float *out=NULL,det;
//  out=dest==NULL?(float*)calloc(9,sizeof(float)):dest;
//  out[0]=a[4]*a[8]-a[7]*a[5];
//  out[1]=a[2]*a[7]-a[1]*a[8];
//  out[2]=a[1]*a[5]-a[2]*a[4];
//  det=out[0]*a[0]+out[1]*a[3]+out[2]*a[6];
//  if(determinant!=NULL) *determinant=det;
//  if(det==0){
//    if(dest==NULL)
//      free(out);
//    return NULL;
//  }
//  out[0]/=det;
//  out[1]/=det;
//  out[2]/=det;
//  out[3]=(a[5]*a[6]-a[3]*a[8])/det;
//  out[4]=(a[0]*a[8]-a[6]*a[2])/det;
//  out[5]=(a[3]*a[2]-a[5]*a[0])/det;
//  out[6]=(a[7]*a[3]-a[6]*a[4])/det;
//  out[7]=(a[6]*a[1]-a[7]*a[0])/det;
//  out[8]=(a[0]*a[4]-a[3]*a[1])/det;
//  return out;
//}
///* computes the product of 3x3 matrix a with vector b*/
//float* mult3x3_3x1(float *a, float *b, float *dest){
//  float c[3],*out;
//  c[0]=a[0]*b[0]+a[3]*b[1]+a[6]*b[2];
//  c[1]=a[1]*b[0]+a[4]*b[1]+a[7]*b[2];
//  c[2]=a[2]*b[0]+a[5]*b[1]+a[8]*b[2];
//  out=dest==NULL?(float*)calloc(3,sizeof(float)):dest;
//  out[0]=c[0];
//  out[1]=c[1];
//  out[2]=c[2];
//  return out;
//}
/*******************************/
/* kernel as defined in eq (20)*/
/*******************************/
float kernelPoly6(float *r, float h){
  static float coef=315./64./PI;
  float roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return 0;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1)
    return 0;
  return coef/(h2*h)*(1-roh2)*(1-roh2)*(1-roh2);
}
float* kernelPoly6Grad(float *r, float h){
  static float coef=-945./32./PI;
  static float grad[3];
  float roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return NULL;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1)
    return NULL;
  grad[0]=grad[1]=grad[2]=coef/(h2*h2*h)*(1-roh2)*(1-roh2);
  grad[0]*=r[0];
  grad[1]*=r[1];
  grad[2]*=r[2];
  return grad;
}
float kernelPoly6Lap(float *r, float h){
  static float coef=-945./32./PI;
  float roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return 0;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1)
    return 0;
  return coef/(h2*h2*h)*(1-roh2)*(3-7*roh2);
}

/*******************************/
/* kernel as defined in eq (21)*/
/*******************************/
float kernelSpiky(float *r, float h){
  static float coef=15./PI;
  float roh,roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return 0;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1)
    return 0;
  roh=sqrt(roh2);
  return coef/(h2*h)*(1-roh)*(1-roh)*(1-roh);
}
float* kernelSpikyGrad(float *r, float h){
  static float coef=-45./PI;
  static float grad[3];
  float roh,roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return NULL;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1 || roh2<0.0001)
    return NULL;
  roh=sqrt(roh2);
  grad[0]=grad[1]=grad[2]=coef/(h2*h2*h)*(roh-2.+1./roh);
  grad[0]*=r[0];
  grad[1]*=r[1];
  grad[2]*=r[2];
  return grad;
}
float kernelSpikyLap(float *r, float h){
  static float coef=-90./PI;
  float roh,roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return 0;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1 || roh2<0.0001)
    return 0;
  roh=sqrt(roh2);
  return coef/(h2*h2*h)*(roh+roh-3.+1./roh);
}

/*******************************/
/* kernel as defined in eq (22)*/
/*******************************/
float kernelVisco(float *r, float h){
  static float coef=15./2./PI;
  float roh,roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return 0;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1. || roh2<0.0001)
    return 0;
  roh=sqrt(roh2);
  return coef/(h2*h)*(-(roh*roh2)/2.+roh2+0.5/roh-1.);
}
float* kernelViscoGrad(float *r, float h){
  static float coef=15./2./PI;
  static float grad[3];
  float roh,roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return NULL;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1. || roh2<0.0001)
    return NULL;
  roh=sqrt(roh2);
  grad[0]=grad[1]=grad[2]=coef/(h2*h2*h)*(-1.5*roh+2.-0.5/(roh*roh2));
  grad[0]*=r[0];
  grad[1]*=r[1];
  grad[2]*=r[2];
  return grad;
}
float kernelViscoLap(float *r, float h){
  static float coef=45./PI;
  float roh2,h2;

  if(fabs(r[0])>=h || fabs(r[1])>=h || fabs(r[2])>=h)
    return 0;
  h2=h*h;
  roh2=(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/h2;
  if(roh2>=1.)
    return 0;
  return coef/(h2*h2*h)*(1.-sqrt(roh2));
}


