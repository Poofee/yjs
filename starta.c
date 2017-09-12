/*
基本初始化元件：
1、得到线性方程组未知量的个数以及每个未知节点对应的方程号；
2、将单元刚度矩阵的每个元素（非零元）与总刚度矩阵元素（非零元）一一对应，对应alch()函数；
3、对总刚度矩阵的每一行进行排序，得到每一行非零元元素的个数与位置，将信
   息存在一维数组na和numcol中，对应blch()函数。
*/
#include "felac.h"
void order(int,int *);
int aclh(int,int *,int *,int,int *,int);
void bclh(int,int *,int *,int *,int *,int *,int);
/* subroutine */
void starta(coor0,dof,nodvar,elem,matrix,f)
int dof,*nodvar;
double **f;
struct coordinates coor0;
struct element elem;
struct matrice *matrix;
{
    /*
    变量解释：
    nelem：  该场某类型单元的单元总数；
    nnode：  某类型单元一个单元上的节点数+1(包括材料)。
    neq：    该场待求方程总数；
    maxcol： 总体刚度矩阵的最大列号；
    maxa:    总刚非零元总数；
    kvar：   总自由度数；
    dim：    坐标维度；
    knode：  节点总数；
    dof：    该场每个节点的自由度数；
    ntype：  该场单元类型数；
    *coor：  节点坐标空间；
    *nodvar：该场节点自由度对应的规格数或者方程号；
    *numcol：总刚矩阵每行(或列)的起始元素在na中的位置；
    *na：    总刚去零元,一维变带宽存储矩阵；
    */
    int nelem,nnode,neq;
    int *numcol,*na,*mht;
    int i,j,l,n,nn,ne,nne,inod,nodi,numel,idof,maxcol,ityp,inv,maxa,kvar;
    int dim,knode,lm[500],node[500];
    double *coor;
    dim   = coor0.dim;
    knode = coor0.knode;
    coor  = coor0.coor;
    kvar  = knode*dof;
    neq=0;
    /*  遍历所有节点上的所有自由度，得到每个未知量对应的方程号以及方程组的方程个数neq  */
    for (j=1; j<=knode; ++j)
        for (i=1; i<=dof; ++i)
            if (nodvar[(i-1)*(knode)+j-1] == 1)
                nodvar[(i-1)*(knode)+j-1] = ++neq;
    /*  从节点未知量的方程号上处理节点的主从约束关系  */
    for (j=1; j<=knode; ++j)
        for (i=1; i<=dof; ++i)
            if (nodvar[(i-1)*(knode)+j-1] < -1)
            {
                n = -nodvar[(i-1)*(knode)+j-1]-1;
                nodvar[(i-1)*(knode)+j-1] = nodvar[(i-1)*(knode)+n-1];
            }
    if (neq>0)
    {
        maxcol=maxt/neq;
    }
    else
    {
        printf("warning! there isn't any unknown variable on nodes!%c",endl);
        maxcol=1;
    }
    na = (int *) calloc(maxt,sizeof(int));
    (*matrix).jdiag=(int *) calloc(neq+1,sizeof(int));
    mht = na;
    numcol = (int *) calloc(neq+1,sizeof(int));
    for (i=0; i<=neq; ++i)
        numcol[i] = 0;
    numel = 0;
    nn = 0;
    /*  遍历每种类型的所有单元，将单元刚度矩阵的每个元素（非零元）与总刚度矩阵的元素（非零元）一一对应  */
    for (ityp=1; ityp<=elem.ntype; ++ityp)
    {
        nelem = elem.nelem[ityp];
        nnode = elem.nnode[ityp];
        for (i=1; i<=nelem; ++i)
        {
            for (j=1; j<=nnode; ++j)
                node[j] = elem.node[++nn];
            nne = nnode-1;
            {
                l = 0;
                /*  在一个单元中，l表示该单元中待求方程数，并将第i方程对应的整体方程号赋给lm[i]  */
                for (inod=1; inod<=nne; ++inod)
                {
                    nodi = node[inod];
                    for (idof=1; idof<=dof; ++idof)
                    {
                        inv = nodvar[(idof-1)*(knode)+nodi-1];
                        if (inv>0)
                        {
                            l++;
                            lm[l] = inv;
                        }
                    }
                }
                numel++;
                /*  aclh函数：得到总刚矩阵非零元的位置(存储在*mht矩阵中)；
                             得到总刚矩阵每行非零元的个数(存储在*numcol中)  */
                if (l>0)
                    aclh(neq,numcol,mht,l,lm,maxcol);
            }
        }
    }
    /*  bclh函数：对mht矩阵按行排序，将其中的非零元素存储在na中；
                  得到总体刚度矩阵每行的起始元素在na中的位置(存储在numcol中)  */
    bclh(neq,numcol,mht,na,(*matrix).jdiag,lm,maxcol);
    /*  把得到的总体刚度矩阵信息存储到matrix中去  */
    maxa=numcol[neq];
    (*matrix).maxa = maxa;
    (*matrix).neq = neq;
    (*matrix).na = (int *)   calloc(maxa+1,sizeof(int));
    (*matrix).a = (double *) calloc(maxa+1,sizeof(double));
    (*matrix).numcol = (int *)  calloc(neq+1,sizeof(int));
    (*f) = (double *) calloc(neq+1,sizeof(double));
    for (n=0; n<=neq; ++n) (*matrix).numcol[n] = numcol[n];
    for (n=0; n<maxa; ++n) (*matrix).na[n+1]   = na[n];
    free(numcol);
    free(na);
    return;
}
