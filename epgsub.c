/* aclh函数：将单元刚度矩阵的每个元素（非零元）与总刚度矩阵元素（非零元）一一对应，存储在mht中
             得到总刚矩阵每行非零元的个数，存储在numcol中 */
/*subroutine*/ int aclh(neq,numcol,mht,nd,lm,maxcol)
int neq,*numcol,*mht,nd,*lm,maxcol;
{
    int i,j,k,ni,nj;
    for (i=1; i<=nd; ++i)
    {
        ni = lm[i];
        for (j=1; j<=nd; ++j)
        {
            nj = lm[j];
            for (k=1; k<=numcol[ni]; ++k)
                if (nj==mht[(ni-1)*maxcol+k-1]) goto l300;
            numcol[ni] = numcol[ni]+1;
            k = numcol[ni];
            mht[(ni-1)*maxcol+k-1] = nj;
            if (numcol[ni]>maxcol)
            {
                // write(*,*) 'exeet maxcol ....',numcol(ni),' > ',maxcol
                return 1;
            }
l300:
            continue;
        }
    }
    return 0;
}
/* bclh函数：对mht矩阵按行排序，将其中的非零元素存储在na中；
             得到总体刚度矩阵每行的起始元素在na中的位置(存储在numcol中) */
/*subroutine*/ void bclh(neq,numcol,mht,na,jdiag,lmi,maxcol)
int neq,*numcol,*mht,*na,*jdiag,*lmi,maxcol;
{
    void order(int,int *);
    int n,nn,i,li,nsum,j;
    for (n=1; n<=neq; ++n)
    {
        li = numcol[n];
        for (i=1; i<=li; ++i)
        {
            lmi[i] = mht[(n-1)*maxcol+i-1];
        }
        order(li,lmi);
        for (i=1; i<=li; ++i)
        {
            mht[(n-1)*maxcol+i-1] = lmi[i];
        }
    }
    nsum = 0;
    for (n=1; n<=neq; ++n)
    {
        for (i=1; i<=numcol[n]; ++i)
        {
            na[nsum] = mht[(n-1)*maxcol+i-1];
            nsum = nsum+1;
        }
    }
    for (n=1; n<=neq; ++n)
    {
        numcol[n] += numcol[n-1];
    }
    for (i=0; i<neq; i++)
    {
        for (j=numcol[i]; j<numcol[i+1]; j++)
        {
            if(na[j]==(i+1)) goto l200;
        }
l200:
        jdiag[i]=j+1;
    }
    return;
}
/* order:排序函数，从小到大排序 */
/*subroutine*/ void order(nd,lm)
int nd,*lm;
{
    int i,j,j0,ls;
    for (i=1; i<=nd; ++i)
    {
        ls = lm[i]+1;
        for (j=i; j<=nd; ++j)
        {
            if (lm[j]<=ls)
            {
                ls = lm[j];
                j0 = j;
            }
        }
        lm[j0] = lm[i];
        lm[i]  = ls;
    }
    return;
}

