#include "particles.h"
#include <string.h>
#ifdef USE_MULTI_THREADING
#include <pthread.h>
#endif
#include <time.h>

//static float getCellKernelCoef(Wall *walls, int nbWalls, float pos[3], float radius);
static int setParticlesGrid(Grid *grid, Particle *part0, int nbParts);

#ifdef USE_MULTI_THREADING
int freeFlightMT(System *sys, float dt, int nbTh);
int checkPPCollisionMT(System *sys, int nbTh);
#endif

/* creates and returns a structure System containing nbParts particles
 * */
System* createSystem(int nbParts, float minMax[3][2]){
  System *sys;
  int i,j;
  float v1[3],v2[3],v3[3];

  if(nbParts<=0) return NULL;
  for(i=0;i<3;i++)
    if(minMax[i][0]>=minMax[i][1])
      return NULL;
  if((sys=(System*)calloc(1,sizeof(System)))==NULL)
    return NULL;
  sys->time=0;
  sys->nbWalls=12;
  sys->walls=(Wall*)calloc(sys->nbWalls,sizeof(Wall));
  for(i=0;i<3;i++){
    for(j=0;j<2;j++)
      sys->boundaries[i][j]=minMax[i][j];
    for(j=0;j<4;j++){
      v1[i]=      minMax[i][(j==0||j==3)?0:1];
      v1[(i+1)%3]=minMax[(i+1)%3][(j==0||j==3)?0:1];
      v1[(i+2)%3]=minMax[(i+2)%3][1-(j>>1)];
      v2[i]=      minMax[i][j&1];
      v2[(i+1)%3]=minMax[(i+1)%3][1-(j&1)];
      v2[(i+2)%3]=minMax[(i+2)%3][1-(j>>1)];
      v3[i]=      minMax[i][(j==0||j==3)?1:0];
      v3[(i+1)%3]=minMax[(i+1)%3][(j==0||j==3)?1:0];
      v3[(i+2)%3]=minMax[(i+2)%3][1-(j>>1)];
      sys->walls[4*i+j]=createWall(v1,v2,v3);
    }
  }
  sys->nbParts=nbParts;
  if((sys->parts=(Particle*)calloc(nbParts,sizeof(Particle)))==NULL)
      return NULL;
  sys->prevStepCollisionTable=(unsigned char*)calloc(nbParts*(nbParts-1)/2,sizeof(unsigned char));
  sys->h=0.4;
  return sys;
}

/* simply free the memory allocated dynamically by sys
 * (sys will also be deleted)*/
int deleteSystem(System *sys){
  free(sys->parts);
  //todo free grid
  free(sys);
  return 0;
}

/*
 * sets the grid used to adjust the value returned by the kernel
 * function when a particle is close to a wall
 */
int setGrid(System *sys, float cellSize){
  Grid *grid;
  Cell *cell, *cell1, **nghb;
//  float cellDim[3];
  int alpha,ixyz[3],index;
  int ix,iy,iz,ixyzNghb[3];
  int gridDim[3];

  grid=(Grid*)calloc(1,sizeof(Grid));
  grid->ncells=1;
  for(alpha=0; alpha<3; alpha++){
    gridDim[alpha]=(int)((sys->boundaries[alpha][1]-sys->boundaries[alpha][0])/cellSize)+1;
    grid->dim[alpha]=gridDim[alpha];
    grid->ncells*=grid->dim[alpha];
    grid->boundaries[alpha][0]=sys->boundaries[alpha][0]-(grid->dim[alpha]*cellSize-(sys->boundaries[alpha][1]-sys->boundaries[alpha][0]))/2;
    grid->boundaries[alpha][1]=grid->boundaries[alpha][0]+grid->dim[alpha]*cellSize;
  }
  grid->cells=(Cell*)calloc(grid->ncells,sizeof(Cell));
  cell1=grid->cells+grid->ncells;
  for(index=0, cell=grid->cells; index<grid->ncells; index++, cell++){
//    printf("%d/%d\n", index, grid->ncells);fflush(stdout);
    ixyz[0]=index%gridDim[0];
    ixyz[1]=(index/gridDim[0])%gridDim[1];
    ixyz[2]=index/(gridDim[0]*gridDim[1]);
    for(alpha=0; alpha<3; alpha++){
      cell->pos[alpha]=cellSize*(ixyz[alpha]+0.5)+grid->boundaries[alpha][0];
      cell->dim[alpha]=cellSize;
    }
    // set the 14 neighbors that will be used to check collisions. (actually only 13 neighbors plus cell itself)
    nghb=cell->nghb;
    for(iz=0; iz<=1; iz++){
      ixyzNghb[2]=ixyz[2]+iz;
      for(iy=(iz==1?-1:0); iy<=1; iy++){
        ixyzNghb[1]=ixyz[1]+iy;
        for(ix=(iz==1?-1:(iy==1?-1:0)); ix<=1; ix++, nghb++){
          ixyzNghb[0]=ixyz[0]+ix;
          for(alpha=0; alpha<3; alpha++)
            if(ixyzNghb[alpha]<0 || ixyzNghb[alpha]>=gridDim[alpha])
              break;
          if(alpha<3)
            *nghb=NULL;
          else
            *nghb=grid->cells+(ixyzNghb[0]+ixyzNghb[1]*gridDim[0]+ixyzNghb[2]*gridDim[0]*gridDim[1]);
        }
      }
    }
    cell->part0=NULL;
  }
  sys->grid=grid;
  printf("setGrid: done\n"); fflush(stdout);
  return 0;
}

/* set the particles in sys.
 * first and last are the indices of the particles which will be set.
 * particles are initialized with zero acceleration and velocity, and
 * randomly distributed within the system.
 * mass is the particle's mass.
 * returns 0 if particles set correctly*/
int setParticles(System *sys, int first, int last, float mass){
  Particle *part,*part0,*part1;
  int      alpha;

  if(first<0 || first>last ||last>=sys->nbParts)
    return 1;
  part0=sys->parts+first;
  part1=sys->parts+last+1;
  for(part=part0;part<part1;part++){
    part->nextPart=NULL;
    part->mass=mass;
    part->radius=0.02;
    part->cor=0.8;
    //part->k=GAS_CONST*TEMP/0.01801528; //for water
    part->k=10; //P=k*ro; P0=1Bar=10^5Pa = k*ro0
    part->k=3; //P=k*ro; P0=1Bar=10^5Pa = k*ro0
    part->ro0=1000.0;                 //for water
    part->mu=1e-3;                    //for water
    part->mu=10;//test
    part->sigma=72.8e-3;              //for water
    for(alpha=0;alpha<3;alpha++){
      part->fg[alpha]=part->fp[alpha]=part->fs[alpha]=part->fv[alpha]=0;
      part->vel[alpha]=0;
      part->pos[alpha]=myrand()*(sys->boundaries[alpha][1]-sys->boundaries[alpha][0]-2.01*part->radius)+sys->boundaries[alpha][0]+1.005*part->radius;
    }
//    part->pos[0]=myrand()*(0.5-2.01*part->radius)+0.5+1.005*part->radius;
//    part->pos[2]=myrand()*(0.5-2.01*part->radius)+0.5+1.005*part->radius;
    part->fg[0]=0;
    part->fg[1]=0;
    part->fg[2]=-9.80665;             //g in m/s^2
  }
  setParticlesGrid(sys->grid, sys->parts, sys->nbParts);
  return 0;
}
/*from the 3 vertices v1, v2, and v3, generate a structure wall*/
Wall createWall(float *v1, float *v2, float *v3){
  Wall  wall;
  float norm,small;
  int   i;

  for(small=i=0;i<3;i++){
    wall.v1[i]=v1[i];
    wall.v2[i]=v2[i];
    wall.v3[i]=v3[i];
    wall.mat[3+i]=v2[i]-v1[i];
    wall.mat[6+i]=v3[i]-v1[i];
    wall.v12[i]=v2[i]-v1[i];//redundant. will have to get rid of mat
    wall.v13[i]=v3[i]-v1[i];//redundant. will have to get rid of mat
    small=small<fabs(wall.v12[i])?fabs(wall.v12[i]):small;
    small=small<fabs(wall.v13[i])?fabs(wall.v13[i]):small;
  }
  small*=1e-5;
  for(norm=i=0;i<3;i++){
    wall.normal[i]=wall.mat[3+(i+1)%3]*wall.mat[6+(i+2)%3]-wall.mat[3+(i+2)%3]*wall.mat[6+(i+1)%3];
    norm+=wall.normal[i]*wall.normal[i];
  }
  norm=sqrt(norm);
  wall.norm=norm;
  wall.dist=0;
  if(norm<small){
    for(i=3;i<9;i++)
      wall.mat[i]=0;
    for(i=0;i<3;i++)
      wall.normal[i]=0;
  }
  else{
    for(i=0;i<3;i++){
      wall.normal[i]/=norm;
      wall.dist+=wall.normal[i]*wall.v1[i];
    }
  }
  return wall;
}
/* add wall in sys
 * */
int addWall(System *sys, Wall wall){
  sys->walls=(Wall*)realloc(sys->walls,(sys->nbWalls+1)*sizeof(Wall));
  sys->walls[sys->nbWalls++]=wall;
  return 0;
}
/* sets the density for each particle, according to equation (3)
 * ASSUMPTION: kernel(r,h)=kernel(-r,h)
 * */
void setDensity(System *sys, float(*kernel)(float*, float)){
  Particle *part1,*part2,*partM;
  float    r[3],h,coef;

  partM=sys->parts+sys->nbParts;
  h=sys->h;
  for(part1=sys->parts;part1<partM;part1++)
    part1->ro=0;
  for(part1=sys->parts;part1<partM;part1++){
    for(part2=part1;part2<partM;part2++){
      r[0]=part1->pos[0]-part2->pos[0];
      r[1]=part1->pos[1]-part2->pos[1];
      r[2]=part1->pos[2]-part2->pos[2];
      coef=kernel(r,h);
      if(coef!=0){
        part1->ro+=part2->mass*coef;
        if(part2!=part1)
          part2->ro+=part1->mass*coef;
      }
    }
  }
//  for(part1=sys->parts;part1<partM;part1++){
//    if(part1->cell && part1->cell->kernelCoef>0)
//      part1->ro/=part1->cell->kernelCoef;
//  }
}
/* sets the force due to the pressure, according to equation (10)
 * ASSUMPTION: kernel(r,h)=-kernel(-r,h)
 * */
void setPressure(System *sys, float*(*kernel)(float*, float)){
  Particle *part1,*part2,*partM;
  float    r[3],h,*coef,coef1,coef2;
  int      alpha;

  partM=sys->parts+sys->nbParts;
  h=sys->h;
  for(part1=sys->parts;part1<partM;part1++){
    part1->p=part1->k*(part1->ro-part1->ro0);
    //part1->p=part1->k*(part1->ro);
    part1->fp[0]=0;
    part1->fp[1]=0;
    part1->fp[2]=0;
  }
  for(part1=sys->parts;part1<partM;part1++){
    for(part2=part1+1;part2<partM;part2++){
      r[0]=part1->pos[0]-part2->pos[0];
      r[1]=part1->pos[1]-part2->pos[1];
      r[2]=part1->pos[2]-part2->pos[2];
      coef=kernel(r,h);
      if(coef!=NULL){
        coef1=part2->mass*(part1->p+part2->p)/(2*part2->ro);
        coef2=part1->mass*(part1->p+part2->p)/(2*part1->ro);
        //coef1=part2->mass*part2->p/part2->ro;
        //coef2=part1->mass*part1->p/part1->ro;
        for(alpha=0;alpha<3;alpha++){
          part1->fp[alpha]-=coef1*coef[alpha];
          if(part2!=part1)
            part2->fp[alpha]+=coef2*coef[alpha]; //+ sign Cf. kernel(-r,h)=-kernel(r,h)
        }
      }
    }
  }
//  for(part1=sys->parts;part1<partM;part1++){
//    if(part1->cell && part1->cell->kernelCoef>0)
//      for(alpha=0;alpha<3;alpha++)
//        part1->fp[alpha]/=part1->cell->kernelCoef;
//  }
}
/* sets the force due to the pressure, according to equation (10)
 * ASSUMPTION: kernel(r,h)=-kernel(-r,h)
 * */
void setPressure2(System *sys, float*(*kernel)(float*, float)){
  Particle *part1,*part2,*partM;
  float    r[3],h,*coef;
  int      alpha;

  partM=sys->parts+sys->nbParts;
  h=sys->h;
  for(part1=sys->parts;part1<partM;part1++){
//    part1->p=part1->k*(part1->ro-part1->ro0);
    part1->p=part1->k*part1->ro;
    part1->fp[0]=0;
    part1->fp[1]=0;
    part1->fp[2]=0;
  }
  for(part1=sys->parts;part1<partM;part1++){
    for(part2=part1+1;part2<partM;part2++){ //kernel(0)=0 (Cf. kernel(-r)=-kernel(r)) so useless to start with part2=part1
      r[0]=part1->pos[0]-part2->pos[0];
      r[1]=part1->pos[1]-part2->pos[1];
      r[2]=part1->pos[2]-part2->pos[2];
      coef=kernel(r,h);
      if(coef!=NULL){
        for(alpha=0;alpha<3;alpha++){
          part1->fp[alpha]-=part2->mass*coef[alpha];
          part2->fp[alpha]+=part1->mass*coef[alpha]; //+ sign Cf. kernel(-r,h)=-kernel(r,h)
        }
      }
    }
  }
  for(part1=sys->parts;part1<partM;part1++){
    for(alpha=0;alpha<3;alpha++)
      part1->fp[alpha]*=part1->k;
//    if(part1->cell && part1->cell->kernelCoef>0)
//      for(alpha=0;alpha<3;alpha++)
//        part1->fp[alpha]/=part1->cell->kernelCoef;
  }
}
/* sets the force due to the viscosity, according to equation (14)
 * ASSUMPTION: kernel(r,h)=kernel(-r,h)
 * */
void setViscosity(System *sys, float(*kernel)(float*, float)){
  Particle *part1,*part2,*partM;
  float    r[3],h,coef,coef1,coef2;
  int      alpha;

  partM=sys->parts+sys->nbParts;
  h=sys->h;
  for(part1=sys->parts;part1<partM;part1++){
    part1->fv[0]=0;
    part1->fv[1]=0;
    part1->fv[2]=0;
  }
  for(part1=sys->parts;part1<partM;part1++){
    for(part2=part1+1;part2<partM;part2++){
      r[0]=part1->pos[0]-part2->pos[0];
      r[1]=part1->pos[1]-part2->pos[1];
      r[2]=part1->pos[2]-part2->pos[2];
      coef=kernel(r,h);
      if(coef!=0){
        coef1=part2->mass*coef/part2->ro;
        coef2=part1->mass*coef/part1->ro;
        for(alpha=0;alpha<3;alpha++){
          part1->fv[alpha]+=coef1*(part2->vel[alpha]-part1->vel[alpha]);
          part2->fv[alpha]+=coef2*(part1->vel[alpha]-part2->vel[alpha]);
        }
      }
    }
    for(alpha=0;alpha<3;alpha++)
      part1->fv[alpha]*=part1->mu;
  }
  for(part1=sys->parts;part1<partM;part1++){
//    if(part1->cell && part1->cell->kernelCoef>0)
//      for(alpha=0;alpha<3;alpha++)
//        part1->fv[alpha]/=part1->cell->kernelCoef;
  }
}
/* sets the force due to the surface tension, according to equation (19)
 * sets also the color field and its gradient
 * ASSUMPTION: kernel(r,h)=kernel(-r,h)
 * */
void setSurfaceT(System *sys, float(*kernel)(float*, float), float*(*kernel2)(float*, float), float(*kernel3)(float*, float)){
  Particle *part1,*part2,*partM;
  float    r[3],h,coef,*coefGrad,coefLap,coef1,coef2,cfLap;
  float    norm,threshold=0.001;
  int      alpha;

  partM=sys->parts+sys->nbParts;
  h=sys->h;
  for(part1=sys->parts;part1<partM;part1++){
    part1->cf=0;
    for(alpha=0;alpha<3;alpha++){
      part1->fs[alpha]=0;
      part1->gradcf[alpha]=0;
    }
  }
  for(part1=sys->parts;part1<partM;part1++){
    for(part2=part1;part2<partM;part2++){
      r[0]=part1->pos[0]-part2->pos[0];
      r[1]=part1->pos[1]-part2->pos[1];
      r[2]=part1->pos[2]-part2->pos[2];
      coef=kernel(r,h);
      if(coef!=0){
        coefGrad=kernel2(r,h);
        coefLap=kernel3(r,h);
        coef1=part2->mass/part2->ro;
        coef2=part1->mass/part1->ro;
        part1->cf+=coef1*coef;
        for(alpha=0;alpha<3;alpha++)
          part1->gradcf[alpha]+=coef1*coefGrad[alpha];
        part1->fs[0]+=coef1*coefLap;
        if(part2!=part1){
          part2->cf+=coef2*coef;
          for(alpha=0;alpha<3;alpha++)
            part2->gradcf[alpha]-=coef2*coefGrad[alpha];
          part2->fs[0]+=coef2*coefLap;
        }
      }
    }
  }
  for(part1=sys->parts;part1<partM;part1++){
    cfLap=part1->fs[0];
    part1->fs[0]=0;
    for(norm=0,alpha=0;alpha<3;alpha++)
      norm+=part1->gradcf[alpha]*part1->gradcf[alpha];
    norm=sqrt(norm);
    if(norm>threshold){
      cfLap*=-part1->sigma/norm;
      for(alpha=0;alpha<3;alpha++)
        part1->fs[alpha]=cfLap*part1->gradcf[alpha];
    }
  }
}
/*simply sets the density and all the forces for each particle
 * */
void setForces(System *sys){
//  setDensity(sys,kernelPoly6);
//  setPressure2(sys,kernelPoly6Grad);
//  setPressure2(sys,kernelSpikyGrad);
  //setPressure(sys,kernelSpikyGrad);
//  setViscosity(sys,kernelViscoLap);
  //setSurfaceT(sys,kernelPoly6,kernelPoly6Grad,kernelPoly6Lap);
}
/*compute forces, at time t, get new velocity based on forces,
 * and update particles position.
 * */
int updateSys(System *sys, float dt){
  float    ft[3];     //total force on particle.
  int      alpha;     //x,y,z index
  Particle *part;
//  float fgs[3][3]={{0,0,-9.81},{9.81,0,0},{0,0,9.81}}; //testnico just for fun
//  float *fg;
//  static float angle=0, gravity;
  static clock_t timeFreeFlight=0;
  static clock_t timeTracking=0;
  static clock_t timePPCollision=0;
  clock_t time1, time2;

  if(sys==NULL){ //testnico
    printf("free flight: %f sec; tracking: %f sec; collision: %f sec\n",
        (float)timeFreeFlight/CLOCKS_PER_SEC, (float)timeTracking/CLOCKS_PER_SEC, (float)timePPCollision/CLOCKS_PER_SEC);
    return 0;
  }

  time1=clock();
#ifdef USE_MULTI_THREADING
  freeFlightMT(sys, dt, THREAD_NUM);
  time2=clock(); timeFreeFlight+=time2-time1; time1=time2;
  setParticlesGrid(sys->grid, sys->parts, sys->nbParts);
  time2=clock(); timeTracking+=time2-time1; time1=time2;
  checkPPCollisionMT(sys, THREAD_NUM);
#else
    for(part=sys->parts; part<sys->parts+sys->nbParts; part++){
      freeFlight(part,dt,sys->walls,sys->nbWalls);
  //    setParticleCell(part, sys);
    }
    time2=clock(); timeFreeFlight+=time2-time1; time1=time2;
    setParticlesGrid(sys->grid, sys->parts, sys->nbParts);
    time2=clock(); timeTracking+=time2-time1; time1=time2;
    checkPPCollision(sys, sys->grid->cells, sys->grid->ncells);
#endif
  time2=clock(); timePPCollision+=time2-time1; time1=time2;

//  angle+=(2*M_PI/4)*dt;
//  if(angle>2*M_PI)
//    angle-=2*M_PI;
//  if(sys->time<1)
//    gravity=0;
//  else
//    gravity=9.81*sin((sys->time-3)*(2*M_PI/4));
//    gravity=9.81;
//    gravity=2;

  for(part=sys->parts; part<sys->parts+sys->nbParts; part++){
    for(alpha=0;alpha<3;alpha++){
      ft[alpha] =part->fp[alpha];
      ft[alpha]+=part->fs[alpha];
      ft[alpha]+=part->fv[alpha];
      //part->vel[alpha]+=dt*(ft[alpha]/part->ro+part->fg[alpha]);
      part->vel[alpha]+=dt*part->fg[alpha];

//      {//testnico just for fun
//        if(sys->time<3)
//          fg=fgs[0];
//        else if(sys->time<6)
//          fg=fgs[1];
//        else
//          fg=fgs[2];
//        part->vel[alpha]+=dt*fg[alpha];
//      }
    }

//    {//testnico just for fun rotating gravity field
//      part->vel[0]+=dt*(part->fg[0]*cos(angle)-part->fg[2]*sin(angle));
//      part->vel[1]+=dt*part->fg[1];
//      part->vel[2]+=dt*(part->fg[0]*sin(angle)+part->fg[2]*cos(angle));
//    }

//    {//testnico just for fun : normal force
//      float norm;
//      for(alpha=0, norm=0; alpha<3; alpha++)
//        norm+=part->pos[alpha]*part->pos[alpha];
//      norm=sqrt(norm);
//      if(norm>0.01){
//        for(alpha=0; alpha<3; alpha++)
//          part->vel[alpha]-=dt*(part->pos[alpha]/norm)*gravity; //+= to repulse, -= to attract
//      }
//    }

//    {//testnico just for fun : tangential force for few seconds
//      if(sys->time<5)
//      {
//        float norm;
//        norm=sqrt(part->pos[0]*part->pos[0]+part->pos[2]*part->pos[2]);
//        if(norm>0.01){
//          part->vel[0]+=dt*(part->pos[2]/norm)*gravity;
//          part->vel[2]-=dt*(part->pos[0]/norm)*gravity;
//        }
//      }
//      else{
//        part->vel[2]-=dt*9.81;
//      }
//    }
  }
  setForces(sys);
  sys->time+=dt;
  sys->prevDt=dt;
  return 0;
}
/*
 * checks if a particle hits a wall during the "free flight"
 * if a collision occurs, updates the position and the velocity.
 * detection of collision is done by computing the intersection
 * of the plane defined by the wall, and the line defined by part->pos
 * and part->vel.
 */
int freeFlight(Particle *part, float dt, Wall *walls, int nbWalls){
  Wall  *wall,*wall1,*wallC;
  float r[3],rr[3],dr[3];
  float cp[3];
  float omega, t=0, beta, gamma, tmin, coef;
  int alpha;

  wall1=walls+nbWalls;
  for(alpha=0; alpha<3; alpha++)
    r[alpha]=part->pos[alpha];
  do{
    wallC=NULL;
    tmin=1;
    for(alpha=0; alpha<3; alpha++)
      dr[alpha]=part->vel[alpha]*dt;
    for(wall=walls;wall<wall1;wall++){
      omega=wall->normal[0]*dr[0]+wall->normal[1]*dr[1]+wall->normal[2]*dr[2];
      if(omega>=0) //particle not going in the good direction
        continue;
      for(alpha=0; alpha<3; alpha++)
        rr[alpha]=r[alpha]-part->radius*wall->normal[alpha]; //rr is the point on the spherical particle that may touch the wall
      t=(wall->normal[0]*(wall->v1[0]-rr[0])+wall->normal[1]*(wall->v1[1]-rr[1])+wall->normal[2]*(wall->v1[2]-rr[2]))/omega;
      if(t>1 || t>=tmin) //t>1 means that the collision will happen in the future
        continue;
      if(t<0){ //t<0 means that the collision happened in the past
        if(t<0.2*part->radius/omega) //due to numerical precision, sometime particles crosses wall. solution: go back in time!
          continue; // here if the particle is less than 10% inside the wall (20% of radius), we go back in time
      }
      // rr+t*dr is on the plane containing wall. have to check now if it is on wall
      omega*=wall->norm;
      cp[0] = (rr[1]-wall->v1[1])*dr[2] - (rr[2]-wall->v1[2])*dr[1];
      cp[1] = (rr[2]-wall->v1[2])*dr[0] - (rr[0]-wall->v1[0])*dr[2];
      cp[2] = (rr[0]-wall->v1[0])*dr[1] - (rr[1]-wall->v1[1])*dr[0];
      beta=-(cp[0]*wall->v13[0]+cp[1]*wall->v13[1]+cp[2]*wall->v13[2])/omega;
      if(beta<-0.001 || beta>1.001)
        continue;
      gamma=(cp[0]*wall->v12[0]+cp[1]*wall->v12[1]+cp[2]*wall->v12[2])/omega;
      if(gamma<-0.001 || gamma>1.001 || (beta+gamma)<-0.001 || (beta+gamma)>1.001)
        continue;
      //collision!
      wallC=wall;
      tmin=t;
    }
    if(wallC!=NULL){//collision detected. update position and velocity
      coef=-(1+part->cor)*(part->vel[0]*wallC->normal[0]+part->vel[1]*wallC->normal[1]+part->vel[2]*wallC->normal[2]);
      for(alpha=0; alpha<3; alpha++){
        r[alpha]+=tmin*dr[alpha]; //update position
        part->vel[alpha]+=coef*wallC->normal[alpha];//reflect velocity
      }
      dt*=1-tmin; //update remaining time
    }
  }while(wallC!=NULL);
  for(alpha=0;alpha<3;alpha++)
    part->pos[alpha]=r[alpha]+dr[alpha];
  return 0;
}

#ifdef USE_MULTI_THREADING
typedef struct _FreeFlightThreadArgs{
  pthread_mutex_t mutex;
  pthread_cond_t condition;
  pthread_mutex_t *masterMutex;
  pthread_cond_t *masterCondition;
  Particle *parts;
  int nbParts;
  float dt;
  Wall *walls;
  int nbWalls;
  int die;
  int done;
} FreeFlightThreadArgs;

void *freeFlightSlave(void *arg){
  FreeFlightThreadArgs *threadArgs;
  Particle *part, *part1;

  threadArgs=(FreeFlightThreadArgs*)arg;
  part1=threadArgs->parts+threadArgs->nbParts;
  pthread_mutex_init(&(threadArgs->mutex), NULL);
  pthread_cond_init(&(threadArgs->condition), NULL);
  pthread_mutex_lock(&(threadArgs->mutex));
  pthread_mutex_lock(threadArgs->masterMutex); //letting master know that thread is ready
  pthread_cond_signal(threadArgs->masterCondition);
  pthread_mutex_unlock(threadArgs->masterMutex);
  while(1){
    pthread_cond_wait(&(threadArgs->condition), &(threadArgs->mutex)); //waiting for a new job
    if(threadArgs->die) //state should be accessed only with threadArgs->mutex locked
      break;
    for(part=threadArgs->parts; part<part1; part++)
      freeFlight(part, threadArgs->dt, threadArgs->walls, threadArgs->nbWalls);
    pthread_mutex_lock(threadArgs->masterMutex); //job done. let the master know
    threadArgs->done=1;
    pthread_cond_signal(threadArgs->masterCondition);
    pthread_mutex_unlock(threadArgs->masterMutex);
  }
  pthread_mutex_unlock(&(threadArgs->mutex));
  pthread_cond_destroy(&(threadArgs->condition));
  pthread_mutex_destroy(&(threadArgs->mutex));
  pthread_mutex_lock(threadArgs->masterMutex);
  pthread_cond_signal(threadArgs->masterCondition);
  pthread_mutex_unlock(threadArgs->masterMutex);
  pthread_exit(0);
  return NULL;
}

int freeFlightMT(System *sys, float dt, int nbTh){
  static int firstTime=1;
  static int nThreads;
  static float deltat;
  static pthread_t *threads=NULL;
  static FreeFlightThreadArgs *threadArgs;
  static pthread_mutex_t masterMutex;
  static pthread_cond_t masterCondition;
  int t, nDone;

  if(sys==NULL){
    if(firstTime==0){
      pthread_mutex_lock(&masterMutex);
      for(t=0; t<nThreads; t++){
        pthread_mutex_lock(&(threadArgs[t].mutex));
        threadArgs[t].die=1;
        pthread_cond_signal(&(threadArgs[t].condition));
        pthread_mutex_unlock(&(threadArgs[t].mutex));
        pthread_cond_wait(&masterCondition, &masterMutex);
      }
      pthread_mutex_unlock(&masterMutex);
      free(threads);
      threads=NULL;
      free(threadArgs);
      threadArgs=NULL;
      nThreads=0;
      firstTime=1;
    }
    return 0;
  }
  if(firstTime){
    int nbPartsSet=0;
    deltat=dt;
    nThreads=nbTh<1?1:nbTh;
    threads=(pthread_t*)calloc(nThreads, sizeof(pthread_t));
    threadArgs=(FreeFlightThreadArgs*)calloc(nThreads, sizeof(FreeFlightThreadArgs));
    pthread_mutex_init(&(masterMutex), NULL);
    pthread_cond_init(&(masterCondition), NULL);
    for(t=0; t<nThreads; t++){
      threadArgs[t].masterCondition=&masterCondition;
      threadArgs[t].masterMutex=&masterMutex;
      threadArgs[t].die=0;
      threadArgs[t].dt=deltat;
      threadArgs[t].walls=sys->walls;
      threadArgs[t].nbWalls=sys->nbWalls;
      threadArgs[t].parts=sys->parts+nbPartsSet;
      if(nbPartsSet<sys->nbParts){
        if(t==nThreads-1)
          threadArgs[t].nbParts=sys->nbParts-nbPartsSet;
        else if(nThreads>sys->nbParts)
          threadArgs[t].nbParts=1;
        else
          threadArgs[t].nbParts=sys->nbParts/nThreads;
      }
      else
        threadArgs[t].nbParts=0; //useless thread...
      nbPartsSet+=threadArgs[t].nbParts;
      pthread_mutex_lock(&masterMutex);
      pthread_create(&(threads[t]), NULL, freeFlightSlave, threadArgs+t);
      pthread_cond_wait(&masterCondition, &masterMutex);
      pthread_mutex_unlock(&masterMutex);
    }
    firstTime=0;
  }
  if(deltat!=dt || (nbTh>=1 && nThreads!=nbTh)){
    freeFlightMT(NULL, 0, 0);
    freeFlightMT(sys, dt, nbTh);
    return 0;
  }
  pthread_mutex_lock(&masterMutex);
  //send signal to start, and wait for all of the threads to be done
  for(t=0; t<nThreads; t++){
    threadArgs[t].done=0;
    pthread_mutex_lock(&(threadArgs[t].mutex));
    pthread_cond_signal(&(threadArgs[t].condition));
    pthread_mutex_unlock(&(threadArgs[t].mutex));
  }
  //wait for threads to be done
  nDone=0;
  while(nDone<nThreads){
    pthread_cond_wait(&masterCondition, &masterMutex);
    for(t=0; t<nThreads; t++){//more than one thread may have sent the signal... I think...
      if(threadArgs[t].done==1){
        nDone++;
        threadArgs[t].done=0;
      }
    }
  }
  pthread_mutex_unlock(&masterMutex);
  return 0;
}
#endif

/*
 * check Particle-Particle collision, and update velocities if needed
 */
int checkPPCollision(System *sys, Cell *cell0, int nCells){
  Particle *part1, *part2;
  float n[3],norm,dotp;
  float coef1,coef2;
  int alpha, inghb;
  int i, i1, i2;
  unsigned char *prevStepCollision;
  Cell *cell, *cell1, *nghbCell;
  float k=500000; //stiffness coef //todo compute k from the size/mass of the particle (100000 looks ok for radius 0.05)

  cell1=cell0+nCells;
  for(cell=cell0; cell<cell1; cell++){
    part1=cell->part0;
    while(part1){
      for(inghb=0; inghb<14; inghb++){
        nghbCell=cell->nghb[inghb];
        if(!nghbCell)
          continue;
        if(nghbCell==cell)
          part2=part1->nextPart;
        else
          part2=nghbCell->part0;
        while(part2){
          for(alpha=norm=0; alpha<3; alpha++){
            n[alpha]=part1->pos[alpha]-part2->pos[alpha];
            norm+=n[alpha]*n[alpha];
          }
          norm=sqrt(norm);
          i1=part1-sys->parts;
          i2=part2-sys->parts;
          if(i1<i2){
            i=i1;
            i1=i2;
            i2=i;
          }
          i=i1*(2*sys->nbParts-3-i1)/2+i2-1;
          prevStepCollision=sys->prevStepCollisionTable+i;
          if(norm<(part1->radius+part2->radius) && norm>0){ //collision detected. norm=0 should never happen...
            if(*prevStepCollision==0){
              for(alpha=dotp=0; alpha<3; alpha++){
                n[alpha]/=norm;
                dotp+=(part2->vel[alpha]-part1->vel[alpha])*n[alpha];
              }
              coef1=dotp*(1+part1->cor)*(part2->mass/(part1->mass+part2->mass));
              coef2=-dotp*(1+part2->cor)*(part1->mass/(part1->mass+part2->mass));
              for(alpha=0; alpha<3; alpha++){
                part1->vel[alpha]+=coef1*n[alpha];
                part2->vel[alpha]+=coef2*n[alpha];
              }
              *prevStepCollision=1;
            }
            else{ //implement some spring force to repulse the two particles, otherwise they form some clusters
              coef1=k*(part1->radius-part1->radius*norm/(part1->radius+part2->radius))*sys->prevDt/part1->mass;
              coef2=-k*(part2->radius-part2->radius*norm/(part1->radius+part2->radius))*sys->prevDt/part1->mass;
              for(alpha=0; alpha<3; alpha++){
                part1->vel[alpha]+=coef1*n[alpha];
                part2->vel[alpha]+=coef2*n[alpha];
              }
            }
          }
          else //no collision. clear flag
            *prevStepCollision=0;
          part2=part2->nextPart;
        }
      }
      part1=part1->nextPart;
    }
  }
  return 0;
}

#ifdef USE_MULTI_THREADING
typedef struct _ThreadArgs{
  pthread_mutex_t mutex;
  pthread_cond_t condition;
  pthread_mutex_t *masterMutex;
  pthread_cond_t *masterCondition;
  int state; //0 <-> idle. 1 <-> working.
  int die;
  int *zStatus;
  System *sys;
  Cell *cell;
  int nCells;
} ThreadArgs;

// thread should send a signal when job is done
// all threads should send the same signal. the master will then have to check which thread is idle.
// the master will send a signal to a specific thread to let him know that he can start working.

void *checkPPCollisionSlave(void *arg){
  ThreadArgs *threadArgs;

  threadArgs=(ThreadArgs*)arg;
  pthread_mutex_init(&(threadArgs->mutex), NULL);
  pthread_cond_init(&(threadArgs->condition), NULL);
  pthread_mutex_lock(&(threadArgs->mutex));
  pthread_mutex_lock(threadArgs->masterMutex); //making sure master waits that thread started
  pthread_cond_signal(threadArgs->masterCondition);//(master must use pthread_cond_wait just after creating the thread)
  pthread_mutex_unlock(threadArgs->masterMutex);
  while(1){
    pthread_cond_wait(&(threadArgs->condition), &(threadArgs->mutex));
    if(threadArgs->die) //state should be accessed only with threadArgs->mutex locked
      break;
    checkPPCollision(threadArgs->sys, threadArgs->cell, threadArgs->nCells);
    //done with this plane. signal master when he is ready
    pthread_mutex_lock(threadArgs->masterMutex);
    threadArgs->state=0;
    *(threadArgs->zStatus)=2;
    pthread_cond_signal(threadArgs->masterCondition);
    pthread_mutex_unlock(threadArgs->masterMutex);
  }
  pthread_mutex_unlock(&(threadArgs->mutex));
  pthread_cond_destroy(&(threadArgs->condition));
  pthread_mutex_destroy(&(threadArgs->mutex));
  pthread_mutex_lock(threadArgs->masterMutex);
  pthread_cond_signal(threadArgs->masterCondition);
  pthread_mutex_unlock(threadArgs->masterMutex);
  pthread_exit(0);
  return NULL;
}

/*
 * check Particle-Particle collision, and update velocities if needed
 * this function uses multithreading. nbTh is the number of working threads.
 */
int checkPPCollisionMT(System *sys, int nbTh){ //MT: Multithreading
  static int firstTime=1;
  static pthread_t  *zThreads=NULL;
  static Cell **zCells=NULL;
  static int nZCells=0;
  static int *zStatus=NULL; //0: plane not done yet. 1: thread is workging on it. 2: plane done
  static ThreadArgs *threadArgs=NULL;
  static int nThreads=0;
  static pthread_mutex_t masterMutex;
  static pthread_cond_t masterCondition;
  int i, t, nZDone;
  int nIdle, iFirstIdle=0;

  if(sys==NULL){
    if(firstTime==0){
      free(zCells);
      zCells=NULL;
      free(zStatus);
      zStatus=NULL;
      pthread_mutex_lock(&masterMutex);
      for(t=0; t<nThreads; t++){
        pthread_mutex_lock(&(threadArgs[t].mutex));
        threadArgs[t].die=1;
        pthread_cond_signal(&(threadArgs[t].condition));
        pthread_mutex_unlock(&(threadArgs[t].mutex));
        pthread_cond_wait(&masterCondition, &masterMutex);
      }
      pthread_mutex_unlock(&masterMutex);
      pthread_mutex_destroy(&masterMutex);
      pthread_cond_destroy(&masterCondition);
      free(threadArgs);
      threadArgs=NULL;
      free(zThreads);
      zThreads=NULL;
      nZCells=0;
      nThreads=0;
      firstTime=1;
    }
    return 0;
  }
  if(firstTime){
    nZCells=sys->grid->dim[2];
    zCells=(Cell**)calloc(nZCells, sizeof(Cell*));
    for(i=0; i<nZCells; i++)
      zCells[i]=sys->grid->cells+i*sys->grid->dim[0]*sys->grid->dim[1];
    zStatus=(int*)calloc(nZCells, sizeof(int));
    pthread_mutex_init(&masterMutex, NULL);
    pthread_cond_init(&masterCondition, NULL);
    nThreads=nbTh<1?1:nbTh;
    zThreads=(pthread_t*)calloc(nThreads, sizeof(pthread_t));
    threadArgs=(ThreadArgs*)calloc(nThreads, sizeof(ThreadArgs));
    for(t=0; t<nThreads; t++){
      threadArgs[t].masterMutex=&masterMutex;
      threadArgs[t].masterCondition=&masterCondition;
      threadArgs[t].state=0;
      threadArgs[t].die=0;
      threadArgs[t].sys=sys;
      threadArgs[t].nCells=sys->grid->dim[0]*sys->grid->dim[1];
      pthread_mutex_lock(&masterMutex);
      pthread_create(&(zThreads[t]), NULL, checkPPCollisionSlave, threadArgs+t);
      pthread_cond_wait(&masterCondition, &masterMutex);
      pthread_mutex_unlock(&masterMutex);
    }
    firstTime=0;
  }
  if(nbTh>=1 && nbTh!=nThreads){ //if number of thread changes, then clean up function, and call it again with the new number of threads.
    checkPPCollisionMT(NULL,0);
    checkPPCollisionMT(sys,nbTh);
    return 0;
  }

  pthread_mutex_lock(&masterMutex);
  for(i=0; i<nZCells; zStatus[i++]=0);
  nZDone=0;
  while(nZDone<nZCells){
    //find an idle thread
    for(t=0, nIdle=0; t<nThreads; t++)
      if(threadArgs[t].state==0)
        if(nIdle++==0)
          iFirstIdle=t;
    if(nIdle==0){ //all threads are busy. wait that one is done and start again
      pthread_cond_wait(&masterCondition, &masterMutex);
      for(t=0; t<nThreads; t++)
        if(threadArgs[t].state==0)
          nZDone++;
      continue;
    }
    //find an XY plane
    for(i=0; i<nZCells; i++){
      if(zStatus[i]==1){
        i++; //a thread is working on this plane. the next plane must not be touched until this thread is done.
        continue;
      }
      if(zStatus[i]==2) //this plane is already done. skip it
        continue;
      if(i>0 && zStatus[i-1]==1) //the plane below and above must not be used by another thread
        continue;
      if(i<(nZCells-1) && zStatus[i+1]==1)
        continue;
      break;
    }
    if(i>=nZCells){ //No plane available. Need to wait that a thread finishes its job and then check again.
      pthread_cond_wait(&masterCondition, &masterMutex);
      for(t=0; t<nThreads; t++)
        if(threadArgs[t].state==0)
          if(nIdle--<=0)
            nZDone++;
      continue;
    }
    //found a plane (#i) and an idle thread (#iFirstIdle).
    pthread_mutex_lock(&(threadArgs[iFirstIdle].mutex));
    threadArgs[iFirstIdle].state=1;
    threadArgs[iFirstIdle].cell=zCells[i];
    threadArgs[iFirstIdle].zStatus=zStatus+i;
    zStatus[i]=1;
    pthread_cond_signal(&(threadArgs[iFirstIdle].condition));
    pthread_mutex_unlock(&(threadArgs[iFirstIdle].mutex));
  }
  pthread_mutex_unlock(&masterMutex);
  return 0;
}

int cleanThreads(){
  freeFlightMT(NULL, 0, 0);
  checkPPCollisionMT(NULL, 0);
  return 0;
}
#endif

/*
 * computes the percentage of the sphere of radius "radius", centered at "pos", which is
 * inside the system defined by "walls"
 * technique used in this function: use a "fine" grid to cut the sphere into small pieces,
 * and check if this piece is on the good side of each wall.
 */
//#define SPHERE_NDIV 21 //number of subdivisions along the diameter. the sphere will fit within a cube of SPHERE_NDIV^3 cells
//static float getCellKernelCoef(Wall *walls, int nbWalls, float pos[3], float radius){
//  float dist;
//  float p[3],vec[3];
//  float sqradius;
//  Wall *wall,*wall1;
//  int nCells;
//  int count,nTotal;
//  float cellDim;
//  float dotpP, dotpPos;
//  float beta, gamma, omega;
//  int i,ix,iy,iz;
//
//  wall1=walls+nbWalls;
//  // check if pos is at a distance greater than radius from any wall first.
//  // in which case coef=1
//  for(wall=walls; wall<wall1; wall++){
//    dist=(pos[0]-wall->v1[0])*wall->normal[0] + (pos[1]-wall->v1[1])*wall->normal[1] + (pos[2]-wall->v1[2])*wall->normal[2];
//    if(dist<radius && dist>0)
//      break;
//  }
//  if(wall==wall1)
//    return 1;
//
//  //there's at least one wall close enough to the sphere
//  nCells=SPHERE_NDIV*SPHERE_NDIV*SPHERE_NDIV;
//  count=nTotal=0;
//  cellDim=2*radius/SPHERE_NDIV;
//  sqradius=radius*radius;
//  for(i=ix=iy=iz=0; i<nCells; i++){
//    if(++ix==SPHERE_NDIV){
//      ix=0;
//      if(++iy==SPHERE_NDIV){
//        iy=0;
//        iz++;
//      }
//    }
//    p[0]=ix*cellDim-radius;
//    p[1]=iy*cellDim-radius;
//    p[2]=iz*cellDim-radius;
//    if(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]>sqradius)
//      continue;
//    nTotal++; //number of cells inside sphere
//    p[0]+=pos[0];
//    p[1]+=pos[1];
//    p[2]+=pos[2];
//    //check if there's a wall between pos and p
//    for(wall=walls; wall<wall1; wall++){
//      dotpP=wall->normal[0]*p[0]+wall->normal[1]*p[1]+wall->normal[2]*p[2];
//      if(dotpP>wall->dist)
//        continue; // final position on wrong side of plane: no collision
//      dotpPos=wall->normal[0]*pos[0]+wall->normal[1]*pos[1]+wall->normal[2]*pos[2];
//      if(dotpPos<wall->dist || dotpPos==dotpP)
//        continue; // initial position on wrong side of plane: no collision
//      vec[0]=(pos[1]-wall->v1[1])*(p[2]-wall->v1[2])-(pos[2]-wall->v1[2])*(p[1]-wall->v1[1]);
//      vec[1]=(pos[2]-wall->v1[2])*(p[0]-wall->v1[0])-(pos[0]-wall->v1[0])*(p[2]-wall->v1[2]);
//      vec[2]=(pos[0]-wall->v1[0])*(p[1]-wall->v1[1])-(pos[1]-wall->v1[1])*(p[0]-wall->v1[0]);
//      omega=(dotpPos-dotpP)*wall->norm;
//      beta=wall->v13[0]*vec[0]+wall->v13[1]*vec[1]+wall->v13[2]*vec[2];
//      if(beta<0 || beta>omega)
//        continue; // ray does not intersect wall
//      gamma=-(wall->v12[0]*vec[0]+wall->v12[1]*vec[1]+wall->v12[2]*vec[2]);
//      if(gamma<0 || gamma>omega || (beta+gamma)>omega)
//        continue; // ray does not intersect wall
//      break; // found a wall between pos and p
//    }
//    if(wall==wall1)
//      count++;
//  }
//  return nTotal?((float)(count)/nTotal):1;
//}

/*
 * Sets the particles in the grid, setting up the linked list in each Cell structure.
 * This function is supposed to be called after updating the particle's position.
 */
static int setParticlesGrid(Grid *grid, Particle *part0, int nbParts){
  float (*boundaries)[2];
  float cellDim[3];
  int   gridDim[3];
  int alpha;
  float *pos;
  int   ixyz[3];
  Particle *part, *part1, **cellPart;
  Cell *initCell, *newCell;

  boundaries=grid->boundaries;
  for(alpha=0; alpha<3; alpha++){
    cellDim[alpha]=grid->cells->dim[alpha];
    gridDim[alpha]=grid->dim[alpha];
  }
  part1=part0+nbParts;
  for(part=part0; part<part1; part++){
    //get index in grid
    pos=part->pos;
    for(alpha=0; alpha<3; alpha++){
      if(pos[alpha]<boundaries[alpha][0] || pos[alpha]>boundaries[alpha][1])
        break;
      ixyz[alpha]=(int)((pos[alpha]-boundaries[alpha][0])/cellDim[alpha]);
      if(ixyz[alpha]<0)               ixyz[alpha]=0;
      if(ixyz[alpha]>=gridDim[alpha]) ixyz[alpha]=gridDim[alpha]-1;
    }
    //check if the particle moved outside of its initial cell
    newCell = alpha<3 ? NULL : grid->cells + (ixyz[0] + (ixyz[1] + ixyz[2]*gridDim[1])*gridDim[0]);
    if(part->cell != newCell){
      //update both initial and final cell
      initCell=part->cell;
      if(initCell){
        cellPart=&(initCell->part0);
        while(*cellPart){
          if(*cellPart==part){
            *cellPart=part->nextPart;
            break;
          }
          cellPart=&((*cellPart)->nextPart);
        }
      }
      if(newCell){
        part->nextPart=newCell->part0;
        newCell->part0=part;
      }
      part->cell=newCell;
    }
  }
  return 0;
}

/*generate a .par file containing the simulation
 * returns 0 if everything went fine.
 *
  the file is built as follow (byte=unsigned char):
  PAR     (3 bytes)
  nbTypes (1 byte)
  name    type 1  (NAME_LENGTH bytes)
  radius  type 1  (1 float)
  red     type 1  (1 byte)
  green   type 1  (1 byte)
  blue    type 1  (1 byte)
  name    type 2  (NAME_LENGTH bytes)
  radius  type 2  (1 float)
  etc.
  nbVariables (1 byte)
  name    x axis     (NAME_LENGTH bytes)
  name    y axis     (NAME_LENGTH bytes)
  name    z axis     (NAME_LENGTH bytes)
  name    variable 1 (NAME_LENGTH bytes)
  name    variable 2 (NAME_LENGTH bytes)
  etc.
  nbParticles (1 unsigned int)
  type of particle 1  (1 byte)
  type of particle 2  (1 byte)
  etc.
  nbBonds       (1 unsigned int)
  bonds radius  (1 float)
  1st index of bond 1 (1 unsigned int)
  2nd index of bond 1 (1 unsigned int)
  1st index of bond 2 (1 unsigned int)
  etc.
  nbSteps       (1 unsigned int)
  time step 1   (1 float)
  x coordinate particle 1 step 1  (1 float)
  y coordinate particle 1 step 1  (1 float)
  z coordinate particle 1 step 1  (1 float)
  variable 1 particle 1 step 1    (1 float)
  variable 2 particle 1 step 1    (1 float)
  etc.
  x coordinate particle 2 step 1  (1 float)
  etc.
  time step 2   (1 float)
  etc.
 * */
int generateParFile(System *sys, char *filename, unsigned char flag){

  static FILE *file=NULL;
  static long posNbSteps=0;
  static long posEnd=0;
  unsigned char tmpc;
  unsigned int  tmpi;
  float         tmpf;
  char          string[NAME_LENGTH];
  int i;

  switch(flag){
  case CREATE:
    if(file!=NULL)
      fclose(file);
    if((file=fopen(filename,"wb+"))==NULL)
      return 1;
    fwrite("PAR",sizeof(char),3,file);
    tmpc=1; //nb types
    fwrite(&tmpc,sizeof(unsigned char),1,file);
    strcpy(string,"type_1");
    fwrite(string,sizeof(char),NAME_LENGTH,file);
    tmpf=sys->parts->radius; //radius
    fwrite(&tmpf,sizeof(float),1,file);
    tmpc=0; //red,green
    fwrite(&tmpc,sizeof(unsigned char),1,file);
    fwrite(&tmpc,sizeof(unsigned char),1,file);
    tmpc=255; //blue
    fwrite(&tmpc,sizeof(unsigned char),1,file);
    tmpc=2; //nb variables
    fwrite(&tmpc,sizeof(unsigned char),1,file);
    strcpy(string,"x [m]");
    fwrite(string,sizeof(char),NAME_LENGTH,file);
    strcpy(string,"y [m]");
    fwrite(string,sizeof(char),NAME_LENGTH,file);
    strcpy(string,"z [m]");
    fwrite(string,sizeof(char),NAME_LENGTH,file);
    strcpy(string,"density [kg/m^3]");
    fwrite(string,sizeof(char),NAME_LENGTH,file);
    strcpy(string,"pressure [Pa]");
    fwrite(string,sizeof(char),NAME_LENGTH,file);
    tmpi=sys->nbParts;
    fwrite(&tmpi,sizeof(unsigned int),1,file);
    tmpc=0; //type of particles
    for(i=0;i<sys->nbParts;i++)
      fwrite(&tmpc,sizeof(unsigned char),1,file);
    tmpi=0; //nb bonds
    fwrite(&tmpi,sizeof(unsigned int),1,file);
    tmpf=0; //bond radius
    fwrite(&tmpf,sizeof(float),1,file);
    posNbSteps=ftell(file);
    tmpi=1; //nb steps
    fwrite(&tmpi,sizeof(unsigned int),1,file);
    tmpf=sys->time; //time first step
    fwrite(&tmpf,sizeof(float),1,file);
    for(i=0;i<sys->nbParts;i++){
      fwrite(sys->parts[i].pos,sizeof(float),3,file);
      fwrite(&(sys->parts[i].ro),sizeof(float),1,file);
      fwrite(&(sys->parts[i].p),sizeof(float),1,file);
    }
    posEnd=ftell(file);
    break;
  case UPDATE:
    if(file==NULL)
      generateParFile(sys,filename,CREATE);
    else{
      fseek(file,posNbSteps,SEEK_SET);
      tmpc=fread(&tmpi,sizeof(unsigned int),1,file);
      tmpi++; //update nb steps
      fseek(file,-sizeof(unsigned int),SEEK_CUR);
      fwrite(&tmpi,sizeof(unsigned int),1,file);
      fseek(file,posEnd,SEEK_SET);
      tmpf=sys->time; //step time
      fwrite(&tmpf,sizeof(float),1,file);
      for(i=0;i<sys->nbParts;i++){
        fwrite(sys->parts[i].pos,sizeof(float),3,file);
        fwrite(&(sys->parts[i].ro),sizeof(float),1,file);
        tmpf =sys->parts[i].vel[0]*sys->parts[i].vel[0]; //compute the norm of the velocity
        tmpf+=sys->parts[i].vel[1]*sys->parts[i].vel[1];
        tmpf+=sys->parts[i].vel[2]*sys->parts[i].vel[2];
        tmpf=sqrt(tmpf);
        fwrite(&tmpf,sizeof(float),1,file);
        //fwrite(&(sys->parts[i].p),sizeof(float),1,file);
      }
      posEnd=ftell(file);
    }
    break;
  case CLOSE:
    if(file!=NULL)
      fclose(file);
    file=NULL;
    posNbSteps=0;
    posEnd=0;
    break;
  default:
    return 1;
  }
  return 0;
}
