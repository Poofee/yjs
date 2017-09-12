/*
������ʼ��Ԫ����
1���õ����Է�����δ֪���ĸ����Լ�ÿ��δ֪�ڵ��Ӧ�ķ��̺ţ�
2������Ԫ�նȾ����ÿ��Ԫ�أ�����Ԫ�����ܸնȾ���Ԫ�أ�����Ԫ��һһ��Ӧ����Ӧalch()������
3�����ܸնȾ����ÿһ�н������򣬵õ�ÿһ�з���ԪԪ�صĸ�����λ�ã�����
   Ϣ����һά����na��numcol�У���Ӧblch()������
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
    �������ͣ�
    nelem��  �ó�ĳ���͵�Ԫ�ĵ�Ԫ������
    nnode��  ĳ���͵�Ԫһ����Ԫ�ϵĽڵ���+1(��������)��
    neq��    �ó����󷽳�������
    maxcol�� ����նȾ��������кţ�
    maxa:    �ܸշ���Ԫ������
    kvar��   �����ɶ�����
    dim��    ����ά�ȣ�
    knode��  �ڵ�������
    dof��    �ó�ÿ���ڵ�����ɶ�����
    ntype��  �ó���Ԫ��������
    *coor��  �ڵ�����ռ䣻
    *nodvar���ó��ڵ����ɶȶ�Ӧ�Ĺ�������߷��̺ţ�
    *numcol���ܸվ���ÿ��(����)����ʼԪ����na�е�λ�ã�
    *na��    �ܸ�ȥ��Ԫ,һά�����洢����
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
    /*  �������нڵ��ϵ��������ɶȣ��õ�ÿ��δ֪����Ӧ�ķ��̺��Լ�������ķ��̸���neq  */
    for (j=1; j<=knode; ++j)
        for (i=1; i<=dof; ++i)
            if (nodvar[(i-1)*(knode)+j-1] == 1)
                nodvar[(i-1)*(knode)+j-1] = ++neq;
    /*  �ӽڵ�δ֪���ķ��̺��ϴ���ڵ������Լ����ϵ  */
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
    /*  ����ÿ�����͵����е�Ԫ������Ԫ�նȾ����ÿ��Ԫ�أ�����Ԫ�����ܸնȾ����Ԫ�أ�����Ԫ��һһ��Ӧ  */
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
                /*  ��һ����Ԫ�У�l��ʾ�õ�Ԫ�д��󷽳�����������i���̶�Ӧ�����巽�̺Ÿ���lm[i]  */
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
                /*  aclh�������õ��ܸվ������Ԫ��λ��(�洢��*mht������)��
                             �õ��ܸվ���ÿ�з���Ԫ�ĸ���(�洢��*numcol��)  */
                if (l>0)
                    aclh(neq,numcol,mht,l,lm,maxcol);
            }
        }
    }
    /*  bclh��������mht���������򣬽����еķ���Ԫ�ش洢��na�У�
                  �õ�����նȾ���ÿ�е���ʼԪ����na�е�λ��(�洢��numcol��)  */
    bclh(neq,numcol,mht,na,(*matrix).jdiag,lm,maxcol);
    /*  �ѵõ�������նȾ�����Ϣ�洢��matrix��ȥ  */
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
