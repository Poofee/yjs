/*
显示最小二乘法算法元件：
利用最小二乘法显式地计算已知量或已求得量的某种表达式的值分布。
*/
#include "felac.h"
void beq4g2(double *,double *,double *,double *,double *,double *,double *,int);
void bet3g2(double *,double *,double *,double *,double *,double *,double *,int);
void beq9g3(double *,double *,double *,double *,double *,double *,double *,int);
void bet6(double *,double *,double *,double *,double *,double *,double *,int);
void ea1b(coor0,dof,elem)
int dof;
struct coordinates coor0;
struct element elem;
{
    /*变量解释：
    dof：    该场每个节点的自由度数；
    ntype：  该场单元类型数；
    nnode:   某类型单元一个单元上的节点数+1(包括材料)；
    kvar:    总自由度数；
    nne：    某类型单元一个单元上的节点数；
    dim：    坐标维数；
    knode：  节点总数；
    neq：    该场待求方程总数；
    nmate：  某种单元类型对应的材料数；
    nprmt：  某种类型单元材料对应的材料参数个数；
    nelem：  该场某类型单元的单元总数；
    mate[]： 材料参数值空间；
    prmt[]:  前nprmt个元素存储单元的材料信息，后nne个元素存储单元的整体节点号；
    coef[]:  耦合场信息；
    node[]： 所有单元的单元拓扑关系空间；
    *nodvar：该场节点自由度对应的规格数或者方程号；
    *coor：  节点坐标空间；
    *r:      单元节点坐标空间
    *es：    单元刚度矩阵；
    *em：    单元质量矩阵；
    *ec：    单元阻尼矩阵；
    *ef：    单元荷载向量。
    */
    int ntype,nnode,kvar;
    int i,j,k,l,m,n,kk,ij,nn,mm,nr,nrw,ne,nne,numel,idof,jdof,
        inod,jnod,nodi,nodj,inv,jnv;
    int neq,node[500],*nodvar;
    double *coor,mate[5000],*r,prmt[500],coef[500];
    double emmax,emmin,*emass;
    int dim,knode;
    int init=0;
    double *es,*em,*ec,*ef;
    double *vectw;
    double *vara1;
    dim   = coor0.dim;
    knode = coor0.knode;
    coor  = coor0.coor;
    kvar  = knode*dof;
    vectw = (double *) calloc(kvar,sizeof(double));
    vara1 = (double *) calloc(knode+1,sizeof(double));
    nodvar = (int *) calloc(kvar,sizeof(int));
    neq = 0;
    /*  遍历所有节点的所有自由度，得到每个未知量对应的方程号以及方程组的方程个数neq  */
    for (j=1; j<=knode; ++j)
        for (i=1; i<=dof; ++i)
            nodvar[(i-1)*(knode)+j-1] = ++neq;
    r = (double *) calloc(500,sizeof(double));
    emass = (double *) calloc(kvar+1,sizeof(double));
    for (n=1; n<=neq; ++n)
        emass[n] = 0.0;
    nrw = 0*dof;
    for (i=1; i<=knode; ++i)
        vara1[i] = unoda[nrw*knode+i-1];
    nrw++ ;
    numel = 0;
    nn = 0;
    mm = 0;
    ntype = nbdetype;
    /* 计算方程的右端项及质量矩阵 */
    for (ityp=1; ityp<=ntype; ++ityp)
    {
        /*  得到单元信息，单元上材料信息，单元上的节点拓扑信息  */
        nmate = elem.nmate[ityp];
        nprmt = elem.nprmt[ityp];
        for (i=1; i<=nmate; ++i)
            for (j=1; j<=nprmt; ++j)
                mate[(i-1)*nprmt+j] = elem.mate[++mm];
        nelem = elem.nelem[ityp];
        nnode = elem.nnode[ityp];
        nne = nnode-1;
        for (ne=1; ne<=nelem; ++ne)
        {
            for (j=1; j<=nnode; ++j)
                node[j] = elem.node[++nn];
            /*  通过每个类型单元的第一个单元的自由度规格数确定该类型单元的单刚，单质等矩阵大小  */
            if (ne==1)
            {
                k=0;
                for (j=1; j<=nne; ++j)
                {
                    jnod = node[j];
                    if (jnod>0)
                        for (l=1; l<=dof; ++l)
                            if (nodvar[(l-1)*(knode)+jnod-1] !=0 )
                                k++;
                }
                kk = k*k;
                es = (double *) calloc(kk,sizeof(double));
                em = (double *) calloc(k+1,sizeof(double));
                ec = (double *) calloc(k+1,sizeof(double));
                ef = (double *) calloc(k+1,sizeof(double));
            }
            /*  读取单元信息：单元的整体节点号(存储在prmt中)，单元节点坐标(存储在r中)  */
            for (j=1; j<=nne; ++j)
            {
                jnod=node[j];
                if (jnod<0) jnod = -jnod;
                prmt[nprmt+j] = jnod;
                i=0;
                coef[j-1+i*nne]=vara1[jnod];
                i++;
                for (i=1; i<=dim; ++i)
                    r[(i-1)*(nne)+j-1] = coor[(i-1)*(knode)+jnod-1];
            }
            imate = node[nnode];
            for (j=1; j<=nprmt; ++j)
                prmt[j] = mate[(imate-1)*nprmt+j];
            switch (ityp)
            {
            case 1 :
                beq4g2(r,coef,prmt,es,em,ec,ef,ne);
                break;
            case 2 :
                bet3g2(r,coef,prmt,es,em,ec,ef,ne);
                break;
            case 3 :
                beq9g3(r,coef,prmt,es,em,ec,ef,ne);
                break;
            case 4 :
                bet6(r,coef,prmt,es,em,ec,ef,ne);
                break;
            }
            i=0;
            for (inod=1; inod<=nne; ++inod)
            {
                nodi=node[inod];
                for (idof=1; idof<=dof; ++idof)
                {
                    inv = nodvar[(idof-1)*(knode)+nodi-1];
                    if (inv==0) goto l600;
                    i++;
                    /* 得到自由度右端项以及质量矩阵 */
                    if (inv>0)
                    {
                        vectw[(idof-1)*(knode)+nodi-1] = vectw[(idof-1)*(knode)+nodi-1]
                                                         +ef[i];
                        emass[inv] = emass[inv]
                                     +em[i];
                    }
l600:
                    continue;
                }
            }
            if (ne==nelem)
            {
                free(es);
                free(em);
                free(ec);
                free(ef);
            }
        }
        numel += nelem;
    }
    free(r);
    /*  保证质量矩阵不为0且不为负数  */
    emmax = 0.0;
    for (i=1; i<=neq; ++i)
        if (emmax<emass[i])
            emmax=emass[i];
    emmin = emmax/1.e008;
    for (i=1; i<=neq; ++i)
        if (fabs(emass[i])<emmin)
            emass[i]=emmin;
    /*  显式地计算未知量  */
    for (i=1; i<=dof; ++i)
        for (j=1; j<=knode; ++j)
        {
            ij = nodvar[(i-1)*(knode)+j-1];
            vectw[(i-1)*(knode)+j-1] /= emass[ij];
        }
//mpi_sere(knode,dof,vectw);
    if (init == 0) init = 1;
// 存储求解的w
    free(unodb);
    n=0;
    n=n+dof;
    unodb = (double *) calloc(n*knode,sizeof(double));
    nrw = 0*dof;
    for (j=1; j<=dof; ++j)
        for (i=1; i<=knode; ++i)
            unodb[(nrw+j-1)*knode+i-1] = vectw[(j-1)*(knode)+i-1];
    nrw += dof;
    free(emass);
    free(vectw);
    free(vara1);
    free(nodvar);
    return;
}
