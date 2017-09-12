#include "felac.h"
int gidpre(void);
int gidmsh(struct coordinates,struct element);
int gidres(struct coordinates);
void starta(struct coordinates,int,int *,struct element, struct matrice *, double **);
void startb(struct coordinates,int,int *,struct element, struct matrice *, double **);
void ea1a(struct coordinates,int,int *,double *,struct element,struct matrice,double *);
void solv(struct matrice,double *);
void ua1a(struct coordinates,int,int *,double *,double *);
void ea1b(struct coordinates,int,struct element);
int main ()
{
    gidpre();
    gidmsh(coor0,elema);
    starta(coor0,dofa,ida,elema,&matrixa,&fa);
    startb(coor0,dofb,idb,elemb,&matrixb,&fb);
    setsolver(default_solver,default_solvpara);
    ea1a(coor0,dofa,ida,ubfa,elema,matrixa,fa);
    solv(matrixa,fa);
    ua1a(coor0,dofa,ida,ubfa,fa);
    setsolver(default_solver,default_solvpara);
    ea1b(coor0,dofb,elemb);
    gidres(coor0);
}
