#ifndef __PARTICLES_H__
#define __PARTICLES_H__

#include <stdlib.h>
#include <stdio.h>
#include "mathe.h"

#define TEMP      (300)            //temperature in Kelvin
#define BOLTZMANN (1.3806504e-23)  //Boltzmann's constant in J/K
#define AVOGADRO  (6.0221418e23)   //Avogadro's constant in 1/mol
#define GAS_CONST (8.314472)       //ideal gas constant in J/K/mol

#define NAME_LENGTH 32
#define CREATE 0
#define UPDATE 1
#define CLOSE  2

#define MAX_PARTS_PER_CELL (16)

#define USE_MULTI_THREADING
#ifdef USE_MULTI_THREADING
#ifndef THREAD_NUM
#define THREAD_NUM (4)
#endif
#endif

typedef struct _Cell{
  float pos[3]; // center coordinates
  float dim[3]; // size along x, y, and z
  struct _Cell* nghb[14];
  struct _Particle* part0;
//  float kernelCoef; // coef used to adjust the kernel
} Cell;

typedef struct _Grid{
  Cell  *cells;
  int   ncells; //=dim[0]*dim[1]*dim[2]
  int   dim[3];
  float boundaries[3][2];
} Grid;

/* structure defining a particle*/
typedef struct _Particle{
  float pos[3];   //position [m]
  float vel[3];   //velocity [m/s]
  float mass;     //mass [kg]
  float radius;
  float cor;      //coefficient or restitution. Used when particle collides with something
  float k;        //gas constant (=GAS_CONST*TEMP/molecular_mass)
  float mu;       //coefficient of viscosity [Pa.s]
  float sigma;    //coefficient of surface tension [N/m]
  float ro0;      //rest density [kg/m^3]
  float ro;       //density [kg/m^3]
  float p;        //pressure
  float cf;       //color field
  float gradcf[3];//gradient of color field
  float fp[3];    //force due to pressure [N]
  float fv[3];    //force due to viscosity [N]
  float fs[3];    //force due to surface tension [N]
  float fg[3];    //force due to gravity [N]
  Cell  *cell;    //cell in which the particle lies
  struct _Particle* nextPart; //used for linked list in Cell structure cell.
} Particle;

/* structure defining a wall=border or obstacle
 * a wall is represented by a triangle, and one
 * side is transparent (particles can go through)
 * while the other side is opaque (particles
 * bounce on it). normal is a vector defining
 * which side is transparent and which one is
 * opaque (same idea as in openGL for example).
 * mat is a 3x3 matrix where the 2nd and 3rd
 * columns are the vectors v1v2 and v1v3
 * respectively.
 * */
typedef struct _Wall{
  float v1[3];      //1st vertex
  float v2[3];      //2nd vertex //can be removed
  float v3[3];      //3rd vertex //can be removed
  float normal[3];  //unit vector=v12xv13/norm
  float dist;       //distance from the origin //can be removed
  float mat[9];     // can be removed

  //vector definition for collision detection
  float v12[3];     //v2-v1
  float v13[3];     //v3-v1
  float norm;       //norm of v12xv13
} Wall;

/* structure defining the system
 * so far contains only an array of particles, but may
 * contain other fields, such as obstacles, etc*/
typedef struct _System{
  float    time;              //time [s]
  float    prevDt;            //previous dt used
  Particle *parts;            //array of particles
  int      nbParts;           //number of particles
  Wall     *walls;            //array of 'walls'
  int      nbWalls;           //number of walls
  float    boundaries[3][2];  //system boundaries ([3] for x/y/z, [2] for min/max)
  float    h;                 //core radius of kernels
  Grid     *grid;
  unsigned char *prevStepCollisionTable; //for testing purpose. used to check previous collision between particles
} System;

System* createSystem(int nbParts, float minMax[3][2]);
int     deleteSystem(System *sys);
int     setParticles(System *sys, int first, int last, float mass);
Wall    createWall(float *v1, float *v2, float *v3);
int     addWall(System *sys,  Wall wall);
int     setGrid(System *sys, float cellSize);
int     freeFlight(Particle *part, float dt, Wall *walls, int nbWalls);
int     checkPPCollision(System *sys, Cell *cell0, int nCells);
void    setDensity(System *sys, float(*kernel)(float*, float));
void    setPressure(System *sys, float*(*kernel)(float*, float));
void    setViscosity(System *sys, float(*kernel)(float*, float));
void    setSurfaceT(System *sys, float(*kernel)(float*, float), float*(*kernel2)(float*, float), float(*kernel3)(float*, float));
void    setForces(System *sys);
int     updateSys(System *sys, float dt);
int     generateParFile(System *sys, char *filename, unsigned char flag);

#ifdef USE_MULTI_THREADING
int     cleanThreads();
#endif

#endif
