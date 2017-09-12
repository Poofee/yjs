#include "felac.h"
#include <string.h>
int getid(FILE *fp,int,int,int *);
int getrd(FILE *fp,int,int,double *);
int getelem(FILE *fp,int,int,struct element *);
int getmate(FILE *fp,struct element *);
int gidpre(void)
{
    FILE *fp;
    int argc,knode;
    char *argv[2];
    int m,k,n,mkd,em,ntype;
    double *coor,d;
    char mystring [160];
    if ((fp = fopen("time0","r"))==NULL)
    {
        printf("cat:cannot open %s\n", "time0");
        /*system("pause");*/ exit(1);
    }
    fgets(mystring,4,fp);
    if(strcmp(mystring,"\xEF\xBB\xBF")!=0)rewind(fp);
    fscanf(fp,"%lf",&tmax);
    fgets(mystring,160,fp);
    fscanf(fp,"%lf",&dt);
    fgets(mystring,160,fp);
    fscanf(fp,"%d",&nstep);
    fgets(mystring,160,fp);
    fscanf(fp,"%d",&itnmax);
    fgets(mystring,160,fp);
    fscanf(fp,"%lf",&tolerance);
    fgets(mystring,160,fp);
    fscanf(fp,"%lf",&dampalfa);
    fgets(mystring,160,fp);
    fscanf(fp,"%lf",&dampbeta);
    fgets(mystring,160,fp);
    fclose(fp);
    itn=1;
    end=0;
    it=0;
    stop=0;
    time_now=0.0;
    setsolver(default_solver,default_solvpara);
    /*
     if ((fp = fopen(argv[1],"r"))==NULL) {
     printf("cat:cannot open %s\n", *argv);
     return 1; }
     fscanf(fp,"%d",&n);
     fgets(mystring,160,fp);
     coor0.knode = n;
     knode = coor0.knode;
    // coor = coor0.coor;
    // getrd(fp,dim,knode,coor);
    */
    if ((fp = fopen("./a1.gid/a1.dat","r"))==NULL)
    {
        printf("cat:cannot open %s\n","./a1.gid/a1.dat");
        /*system("pause");*/ exit(1);
    }
    fscanf(fp,"%d",&n);
    fgets(mystring,160,fp);
    coor0.knode = n;
    knode = n;
    dim = 2;
    coor0.dim = dim;
    k = knode*dim;
    coor0.coor = (double *) calloc(k,sizeof(double));
    coor = coor0.coor;
    getrd(fp,dim,knode,coor);
    nbdetype = 4;
    dofa = 1;
    inita = 0;
    mkd = knode*dofa;
    ubfa = (double *) calloc(mkd,sizeof(double));
    ida = (int *) calloc(mkd,sizeof(int));
    unoda = (double *) calloc(mkd*3,sizeof(double));
    for (m=0; m<(mkd*3); m++) unoda[m] = 0.0;
    dofb = 2;
    initb = 0;
    mkd = knode*dofb;
    ubfb = (double *) calloc(mkd,sizeof(double));
    idb = (int *) calloc(mkd,sizeof(int));
    unodb = (double *) calloc(mkd*3,sizeof(double));
    for (m=0; m<(mkd*3); m++) unodb[m] = 0.0;
    m = getid(fp,dofa,knode,ida);
    m = getrd(fp,dofa,knode,ubfa);
    m = getid(fp,dofb,knode,idb);
    m = getrd(fp,dofb,knode,ubfb);
    fgets(mystring,160,fp);
    fscanf(fp,"%d",&m);
    em = 1-m;
    ntype = 4;
    m = getelem(fp,ntype,em,&elema);
    em = 1- m;
    ntype = 4;
    m = getelem(fp,ntype,em,&elemb);
    em = 1- m;
    fclose(fp);
    if ((fp = fopen("a1.mat","r"))==NULL)
    {
        printf("cat:cannot open %s\n","a1.mat");
        /*system("pause");*/ exit(1);
    }
    m = getmate(fp,&elema);
    m = getmate(fp,&elemb);
    fclose(fp);

    return 0;
}
int getrd(fp,n,m,rp)
FILE *fp;
int n,m;
double *rp;
{
    int i,j,k;
    double d;
    char mystring [160];
    for (j=0; j<m*n; j++) rp[j] = 0.0;
    k = 1;
    fgets(mystring,160,fp);
    while (k>0)
    {
        fscanf(fp,"%d",&k);
        if (k>0)
        {
            for (i=0; i<n; i++)
            {
                fscanf(fp,"%lf",&d);
                rp[m*i+k-1] = d;
            }
        }
    }
    return k;
}

int getid(fp,n,m,ip)
FILE *fp;
int n,m;
int *ip;
{
    int i,j,k;
    int d;
    char mystring [160];
    for (j=0; j<m*n; j++) ip[j] = 1;
    k = 1;
    fgets(mystring,160,fp);
    while (k>0)
    {
        fscanf(fp,"%d",&k);
        if (k>0)
        {
            for (i=0; i<n; i++)
            {
                fscanf(fp,"%d",&d);
                ip[m*i+k-1] = d;
            }
        }
    }
    return k;
}

int getelem(FILE *fp,int ntype,int em,struct element *elem)
{
    int i,j,k,m,n;
    int d,*node;
    int l;
    char mystring [160];
    node = (int *) calloc(maxnode,sizeof(int));
    (*elem).ntype = ntype;
    m=em;
    (*elem).nnode[1]=m;
//  fgets(mystring,160,fp);
    n=0;
    l=0;
    j=1;
    while (j<=ntype)
    {
        fgets(mystring,160,fp);
        fscanf(fp,"%d",&k);
        if (k<0)
        {
            m = 1-k;
            j++;
            (*elem).nnode[j]=m;
            l=0;
        }
        else
        {
            for (i=0; i<m; i++)
            {
                fscanf(fp,"%d",&d);
                n++;
                node[n-1] = d;
            }
            l++;
            (*elem).nelem[j]=l;
        }
    }
    (*elem).node = (int *) calloc(n+1,sizeof(int));
    for (j=0; j<n; ++j)
    {
        (*elem).node[j+1] = node[j];
    }
    free(node);
    return k;
}

int getmate(FILE *fp,struct element *elem)
{
    int n,i,j,m,k;
    double r,*mate;
    char mystring [160];
    mate = (double *) calloc(maxmate,sizeof(double));
    m = 0;
    fscanf(fp,"%d",&k);
    (*elem).ntype = k;
    fgets(mystring,160,fp);
    for (n=1; n<=(*elem).ntype; n++)
    {
        fscanf(fp,"%d",&k);
        (*elem).nmate[n] = k;
        fscanf(fp,"%d",&k);
        (*elem).nprmt[n] = k;
        for (i=1; i<=(*elem).nmate[n]; i++)
        {
            for (j=1; j<=(*elem).nprmt[n]; j++)
            {
                fscanf(fp,"%lf",&r);
                mate[++m] = r;
            }
//           do c=getc(fp); while (c|='\n');
            fgets(mystring,160,fp);
        }
    }
    (*elem).mate = (double *) calloc(m+1,sizeof(double));
    for (i=1; i<=m; ++i) (*elem).mate[i] = mate[i];
    free(mate);
    return m;
}

void chid(coor0,dof,nodvar,elem)
struct coordinates coor0;
int dof,*nodvar;
struct element elem;
{
    int numel,nn,ityp,nelem,nnode,nne,n0;
    int i,j,node[500],inod,nodi,kvar;
    int dim,knode;
    double *coor;
    dim = coor0.dim;
    knode = coor0.knode;
    coor = coor0.coor;
    kvar = knode*dof;
    numel=0;
    nn=0;
    for (ityp=1; ityp<=elem.ntype; ++ityp)
    {
        /*2000*/
        nelem=elem.nelem[ityp];
        nnode=elem.nnode[ityp];
        nne = nnode-1;
        n0=nne;
        if (nne==6) n0=3;
        if (nne==9) n0=4;
        if (nne==10) n0=4;
        if (nne==27) n0=8;
        if (n0!=nne)
        {
            for (i=1; i<=nelem; ++i)
            {
                for (j=1; j<=nnode; ++j) node[j]=elem.node[++nn];
                for (inod=n0+1; inod<=nne; ++inod)
                {
                    nodi=node[inod];
                    nodvar[(dof-1)*knode+nodi-1]=0;
                }
            }
        }
    } /*2000*/
}

