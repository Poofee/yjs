#include "felac.h"
int gidmsh(struct coordinates,struct element);
int gidres(coor0)
struct coordinates coor0;
{
    FILE *fp;
    int i,j,dim,knode,ik;
// double *coor,d;
    int *mlgnode;
    char str[100];
    knode = coor0.knode;
    dim = coor0.dim;
    if (it==0) it=1;
    if ((it%nstep)!=0) goto l100;
    ik = it/nstep;
    if ((fp = fopen("./a1.gid/a1.post.res","a"))==NULL)
    {
        printf("cat:cannot open %s\n", "a1.post.res");
        return  1;
    }
    fprintf(fp,"%s%d%s\n","Result \"unoda0\" \"Load Analysis\"  ",ik," Scalar OnNodes");
    fprintf(fp,"%s\n","ComponentNames \"az\" ");
    fprintf(fp,"%s\n","Values");
    for (j=1; j<=knode; ++j)
    {
        fprintf(fp,"%d",j );
        fprintf(fp," %e",unoda[0*knode+j-1]);
        fprintf(fp,"\n");
    }
    fprintf(fp,"%s\n","End Values");
    fprintf(fp,"%s%d%s\n","Result \"unodb0\" \"Load Analysis\"  ",ik," Vector OnNodes");
    fprintf(fp,"%s\n","ComponentNames \"hx\" \"hy\" ");
    fprintf(fp,"%s\n","Values");
    for (j=1; j<=knode; ++j)
    {
        fprintf(fp,"%d",j );
        fprintf(fp," %e",unodb[0*knode+j-1]);
        fprintf(fp," %e",unodb[1*knode+j-1]);
        fprintf(fp,"\n");
    }
    fprintf(fp,"%s\n","End Values");
    fclose(fp);
l100:
    return 0;
}
