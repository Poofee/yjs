/*
��ʾ��С���˷��㷨Ԫ����
������С���˷���ʽ�ؼ�����֪�������������ĳ�ֱ��ʽ��ֵ�ֲ���
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
    /*�������ͣ�
    dof��    �ó�ÿ���ڵ�����ɶ�����
    ntype��  �ó���Ԫ��������
    nnode:   ĳ���͵�Ԫһ����Ԫ�ϵĽڵ���+1(��������)��
    kvar:    �����ɶ�����
    nne��    ĳ���͵�Ԫһ����Ԫ�ϵĽڵ�����
    dim��    ����ά����
    knode��  �ڵ�������
    neq��    �ó����󷽳�������
    nmate��  ĳ�ֵ�Ԫ���Ͷ�Ӧ�Ĳ�������
    nprmt��  ĳ�����͵�Ԫ���϶�Ӧ�Ĳ��ϲ���������
    nelem��  �ó�ĳ���͵�Ԫ�ĵ�Ԫ������
    mate[]�� ���ϲ���ֵ�ռ䣻
    prmt[]:  ǰnprmt��Ԫ�ش洢��Ԫ�Ĳ�����Ϣ����nne��Ԫ�ش洢��Ԫ������ڵ�ţ�
    coef[]:  ��ϳ���Ϣ��
    node[]�� ���е�Ԫ�ĵ�Ԫ���˹�ϵ�ռ䣻
    *nodvar���ó��ڵ����ɶȶ�Ӧ�Ĺ�������߷��̺ţ�
    *coor��  �ڵ�����ռ䣻
    *r:      ��Ԫ�ڵ�����ռ�
    *es��    ��Ԫ�նȾ���
    *em��    ��Ԫ��������
    *ec��    ��Ԫ�������
    *ef��    ��Ԫ����������
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
    /*  �������нڵ���������ɶȣ��õ�ÿ��δ֪����Ӧ�ķ��̺��Լ�������ķ��̸���neq  */
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
    /* ���㷽�̵��Ҷ���������� */
    for (ityp=1; ityp<=ntype; ++ityp)
    {
        /*  �õ���Ԫ��Ϣ����Ԫ�ϲ�����Ϣ����Ԫ�ϵĽڵ�������Ϣ  */
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
            /*  ͨ��ÿ�����͵�Ԫ�ĵ�һ����Ԫ�����ɶȹ����ȷ�������͵�Ԫ�ĵ��գ����ʵȾ����С  */
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
            /*  ��ȡ��Ԫ��Ϣ����Ԫ������ڵ��(�洢��prmt��)����Ԫ�ڵ�����(�洢��r��)  */
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
                    /* �õ����ɶ��Ҷ����Լ��������� */
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
    /*  ��֤��������Ϊ0�Ҳ�Ϊ����  */
    emmax = 0.0;
    for (i=1; i<=neq; ++i)
        if (emmax<emass[i])
            emmax=emass[i];
    emmin = emmax/1.e008;
    for (i=1; i<=neq; ++i)
        if (fabs(emass[i])<emmin)
            emass[i]=emmin;
    /*  ��ʽ�ؼ���δ֪��  */
    for (i=1; i<=dof; ++i)
        for (j=1; j<=knode; ++j)
        {
            ij = nodvar[(i-1)*(knode)+j-1];
            vectw[(i-1)*(knode)+j-1] /= emass[ij];
        }
//mpi_sere(knode,dof,vectw);
    if (init == 0) init = 1;
// �洢����w
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
