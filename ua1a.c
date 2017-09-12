/*
��������������Ԫ����
1�����������䣬��ʩ����ǿ���������ı߽��ϵ����ɶ����䵽�����������Ľ������У�
2�����Ƿ��������⣬��Խ���е����������ǰ�������Ľ����
3���ж��Ƿ��������������������ʶ��end��1��
*/
#include "felac.h"
/* subroutine */
void ua1a(coor0,dof,nodvar,u,f)
int dof,*nodvar;
double *u,*f;
struct coordinates coor0;
{
    /*
    �������ͣ�
     dim��   ����ά��
     knode�� �ڵ�����
     kvar:   �����ɶ�����
     dof��   �ó�ÿ���ڵ�����ɶ�����
    *coor��  �ڵ�����ռ䣻
    *nodvar���ó��ڵ����ɶȶ�Ӧ�Ĺ�������߷��̺ţ�
    */
    int i,j,n,nrw,kvar;
    int dim,knode;
    double *coor;
    double *vectu;
    dim   = coor0.dim;
    knode = coor0.knode;
    coor  = coor0.coor;
    kvar  = knode*dof;
    vectu = (double *) calloc(kvar,sizeof(double));
    /*  ���������䣬��ʩ����ǿ���������ı߽��ϵ����ɶ����䵽�����������Ľ������У� */
    for (i=1; i<=dof; ++i)
        for (j=1; j<=knode; ++j)
        {
            vectu[(i-1)*(knode)+j-1] = u[(i-1)*(knode)+j-1];
            n = nodvar[(i-1)*(knode)+j-1];
            if (n>0)
                vectu[(i-1)*(knode)+j-1] = f[n];
        }
//mpi_sere(knode,dof,vectu);
// �洢����u
    free(unoda);
    n=0;
    n=n+dof;
    unoda = (double *) calloc(n*knode,sizeof(double));
    nrw = 0*dof;
    for (j=1; j<=dof; ++j)
        for (i=1; i<=knode; ++i)
            unoda[(nrw+j-1)*knode+i-1] = vectu[(j-1)*(knode)+i-1];
    nrw += dof;
    free(vectu);
    return;
}