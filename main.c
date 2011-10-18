/*
 * This code tries to implement a particle-based fluid simulation
 * based on the paper "Particle-Based Fluid Simulation for Interactive Applications"
 * form M. Muller, D. Charypar and M. Gross (SIGGRAPH'05 Symposium on Computer Animation (SCA 2005)).
 * There is no interaction or graphic interface implemented. The program generate a .par file
 * containing the evolution of the particles through the simulation. (a .par file is supposed to
 * be open with the Basic Particle Viewer program (bpv))
 * */

#include "particles.h"
#include <time.h>

void tmpAddWalls_1(System *sys);
void tmpAddWalls_2(System *sys);
void tmpAddWalls_3(System *sys);
//void printGrid(System *sys);

int main(int argc, char *argv[]){

  int    nbParts=20000;
  float  minMax[3][2]={{-1,1},{-1,1},{-1,1}};
  float  cellSize=2*0.02*1.001;// 2*particle radius + 0.1%
  float  dt=0.002,simTime=5;
  float  dtFrame=0.02; //0.04 <=> 25 frames per seconds
  int    nbDtPerFrame,i;
  char   filename[]="toto.par";
  float  mass;
  System *sys;
  clock_t time;

/*
  FILE *file;
  float r[3];

  file=fopen("filter.plt","w");
  fprintf(file,"VARIABLES = x w\n");
  r[1]=r[2]=0;
  for(i=0;i<=100;i++){
    r[0]=i/100.;
    fprintf(file,"%e\t%e\n",i/100.,kernelPoly6(r,1));
  }
  fclose(file);
  return 0;*/


  nbDtPerFrame=dtFrame/dt;
  sys=createSystem(nbParts,minMax);
  tmpAddWalls_1(sys);
  setGrid(sys, cellSize);
  //printGrid(sys); //test nico
  //return 0;
  mass=1000*(minMax[0][1]-minMax[0][0])*(minMax[1][1]-minMax[1][0])*(minMax[2][1]-minMax[2][0])/nbParts;
  setParticles(sys,0,nbParts-1,mass/4); //for water, ro0=1000kg/m3, so 1000 parts/Liter => mi=1e-3kg
  setForces(sys);
  generateParFile(sys,filename,CREATE);

  i=0;
  time=clock();
  while(sys->time<=simTime){
    printf("time = %e (i=%d)\n",sys->time,i);fflush(stdout);
    updateSys(sys,dt);
    if((++i)%nbDtPerFrame==0)
      generateParFile(sys,filename,UPDATE);
  }
  time=clock()-time;

  generateParFile(sys,filename,CLOSE);
  deleteSystem(sys);

  updateSys(NULL,0);
#ifdef USE_MULTI_THREADING
  cleanThreads();
#endif
  printf("done in %f sec!\n", (double)(time)/CLOCKS_PER_SEC);
  return 0;
}

/******************************************************************************************
 * ****************************************************************************************
 * ****************************************************************************************/

void tmpAddWalls_1(System *sys){
  float v1[3],v2[3],v3[3],v4[3];
  v1[0]=-0.5;  v1[1]=1;   v1[2]=-0.5;
  v2[0]=-0.5;  v2[1]=-1;  v2[2]=-0.5;
  v3[0]=1;     v3[1]=-1;  v3[2]=0.5;
  v4[0]=1;     v4[1]=1;   v4[2]=0.5;
  addWall(sys,createWall(v1,v2,v3));
  addWall(sys,createWall(v1,v3,v4));

  v1[0]=-0.5;  v1[1]=1;   v1[2]=-0.51;
  v2[0]=-0.5;  v2[1]=-1;  v2[2]=-0.51;
  v3[0]=1;     v3[1]=-1;  v3[2]=0.49;
  v4[0]=1;     v4[1]=1;   v4[2]=0.49;
  addWall(sys,createWall(v1,v3,v2));
  addWall(sys,createWall(v1,v4,v3));
}

void tmpAddWalls_2(System *sys){
  float minx,maxx;
  float miny,maxy;
  float minz,maxz;
  float v1[3],v2[3],v3[3],v4[3];
  float v5[3],v6[3],v7[3],v8[3];
  float v9[3];

  minx=-1;
  maxx=1;
  miny=-1;
  maxy=1;
  minz=-1;
  maxz=0;

  v1[0]=(minx+maxx)/2;  v1[1]=miny;           v1[2]=(minz+maxz)/2;
  v2[0]=maxx;           v2[1]=(miny+maxy)/2;  v2[2]=(minz+maxz)/2;
  v3[0]=(minx+maxx)/2;  v3[1]=maxy;           v3[2]=(minz+maxz)/2;
  v4[0]=minx;           v4[1]=(miny+maxy)/2;  v4[2]=(minz+maxz)/2;

  v5[0]=minx;  v5[1]=miny;  v5[2]=maxz;
  v6[0]=maxx;  v6[1]=miny;  v6[2]=maxz;
  v7[0]=maxx;  v7[1]=maxy;  v7[2]=maxz;
  v8[0]=minx;  v8[1]=maxy;  v8[2]=maxz;

  v9[0]=(minx+maxx)/2;  v9[1]=(miny+maxy)/2;  v9[2]=minz;

  addWall(sys,createWall(v1,v2,v9));
  addWall(sys,createWall(v1,v6,v2));
  addWall(sys,createWall(v2,v3,v9));
  addWall(sys,createWall(v2,v7,v3));
  addWall(sys,createWall(v3,v4,v9));
  addWall(sys,createWall(v3,v8,v4));
  addWall(sys,createWall(v4,v1,v9));
  addWall(sys,createWall(v4,v5,v1));
}

void tmpAddWalls_3(System *sys){
  float xmin,xmax;
  float ymin,ymax;
  float zmin,zmax; //position of wall
  float holeDim; //hole centered on (xmax+xmin)/2 and (ymin+ymax)/2
  float xmin2,xmax2;
  float ymin2,ymax2;
  float v[24][3];

  xmin=-1;
  xmax=1;
  ymin=-1;
  ymax=1;
  zmin=-0.3;
  zmax=0.3;
  holeDim=0.2;

  xmin2=(xmin+xmax)/2-holeDim/2;
  xmax2=(xmin+xmax)/2+holeDim/2;
  ymin2=(ymin+ymax)/2-holeDim/2;
  ymax2=(ymin+ymax)/2+holeDim/2;

  v[0][0]=xmin;    v[0][1]=ymin;   v[0][2]=zmin;
  v[1][0]=xmin;    v[1][1]=ymax;   v[1][2]=zmin;
  v[2][0]=xmax;    v[2][1]=ymax;   v[2][2]=zmin;
  v[3][0]=xmax;    v[3][1]=ymin;   v[3][2]=zmin;
  v[4][0]=xmin2;   v[4][1]=ymin2;  v[4][2]=zmin;
  v[5][0]=xmin2;   v[5][1]=ymax2;  v[5][2]=zmin;
  v[6][0]=xmax2;   v[6][1]=ymax2;  v[6][2]=zmin;
  v[7][0]=xmax2;   v[7][1]=ymin2;  v[7][2]=zmin;
  v[8][0]=xmin2;   v[8][1]=ymin;   v[8][2]=zmin;
  v[9][0]=xmin2;   v[9][1]=ymax;   v[9][2]=zmin;
  v[10][0]=xmax2;  v[10][1]=ymax;  v[10][2]=zmin;
  v[11][0]=xmax2;  v[11][1]=ymin;  v[11][2]=zmin;
  v[12][0]=xmin;   v[12][1]=ymin;  v[12][2]=zmax;
  v[13][0]=xmin;   v[13][1]=ymax;  v[13][2]=zmax;
  v[14][0]=xmax;   v[14][1]=ymax;  v[14][2]=zmax;
  v[15][0]=xmax;   v[15][1]=ymin;  v[15][2]=zmax;
  v[16][0]=xmin2;  v[16][1]=ymin2; v[16][2]=zmax;
  v[17][0]=xmin2;  v[17][1]=ymax2; v[17][2]=zmax;
  v[18][0]=xmax2;  v[18][1]=ymax2; v[18][2]=zmax;
  v[19][0]=xmax2;  v[19][1]=ymin2; v[19][2]=zmax;
  v[20][0]=xmin2;  v[20][1]=ymin;  v[20][2]=zmax;
  v[21][0]=xmin2;  v[21][1]=ymax;  v[21][2]=zmax;
  v[22][0]=xmax2;  v[22][1]=ymax;  v[22][2]=zmax;
  v[23][0]=xmax2;  v[23][1]=ymin;  v[23][2]=zmax;

  addWall(sys,createWall(v[0],v[1],v[9]));
  addWall(sys,createWall(v[0],v[9],v[8]));
  addWall(sys,createWall(v[5],v[9],v[10]));
  addWall(sys,createWall(v[5],v[10],v[6]));
  addWall(sys,createWall(v[8],v[4],v[7]));
  addWall(sys,createWall(v[8],v[7],v[11]));
  addWall(sys,createWall(v[11],v[10],v[2]));
  addWall(sys,createWall(v[11],v[2],v[3]));

  addWall(sys,createWall(v[12],v[21],v[13]));
  addWall(sys,createWall(v[12],v[20],v[21]));
  addWall(sys,createWall(v[17],v[22],v[21]));
  addWall(sys,createWall(v[17],v[18],v[22]));
  addWall(sys,createWall(v[20],v[19],v[16]));
  addWall(sys,createWall(v[20],v[23],v[19]));
  addWall(sys,createWall(v[23],v[14],v[22]));
  addWall(sys,createWall(v[23],v[15],v[14]));

  addWall(sys,createWall(v[4],v[16],v[19]));
  addWall(sys,createWall(v[4],v[19],v[7]));
  addWall(sys,createWall(v[5],v[18],v[17]));
  addWall(sys,createWall(v[5],v[6],v[18]));
  addWall(sys,createWall(v[5],v[17],v[16]));
  addWall(sys,createWall(v[5],v[16],v[4]));
  addWall(sys,createWall(v[6],v[19],v[18]));
  addWall(sys,createWall(v[6],v[7],v[19]));
}

//void printGrid(System *sys){
//  char filename[]="grid.par";
//  System sysGrid;
//  int i;
//
//  sysGrid.time=0;
//  sysGrid.nbParts=sys->grid->ncells;
//  sysGrid.parts=(Particle*)calloc(sysGrid.nbParts, sizeof(Particle));
//  for(i=0; i<sysGrid.nbParts; i++){
//    sysGrid.parts[i].pos[0]=sys->grid->cells[i].pos[0];
//    sysGrid.parts[i].pos[1]=sys->grid->cells[i].pos[1];
//    sysGrid.parts[i].pos[2]=sys->grid->cells[i].pos[2];
//    sysGrid.parts[i].ro=sys->grid->cells[i].kernelCoef;
//  }
//  generateParFile(&sysGrid,filename,CREATE);
//  generateParFile(&sysGrid,filename,CLOSE);
//  free(sysGrid.parts);
//  printf("generated %s\n",filename);
//  fflush(stdout);
//}
