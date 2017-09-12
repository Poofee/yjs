#include "felac.h"
int gidmsh(coor0,elem)
struct coordinates coor0;
struct element elem;
{
    FILE *fp;
    int i,j,n,dim,knode,*nnode,*node,*nelem;
    double *coor,d;
    int *mlgnode;
    char str[100];
    knode = coor0.knode;
    dim = coor0.dim;
    coor = coor0.coor;
    nnode = elem.nnode;
    node = elem.node;
    nelem = elem.nelem;
    n = 0;
    if ((fp = fopen("./a1.gid/a1.post.msh","w"))==NULL)
    {
        printf("cat:cannot open %s\n", "a1.post.msh");
        return  1;
    }
    fprintf(fp,"%s\n","Mesh \"aeq4g2\" Dimension 2 Elemtype Quadrilateral Nnode 4");
    fprintf(fp,"%s\n","Coordinates");
    for (j=1; j<=knode; ++j)
    {
        fprintf(fp,"%d",j);
        for (i=1; i<=dim; ++i)
        {
            fprintf(fp," %e",coor[(i-1)*knode+j-1]);
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"%s\n","End coordinates");
    fprintf(fp,"%s\n","Elements");
    for (i=1; i<=nelem[1]; ++i)
    {
        fprintf(fp,"%d",i);
        for (j=1; j<=nnode[1]; ++j)
        {
            fprintf(fp," %d",node[++n]);
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"%s\n","End elements");
    fprintf(fp,"%s\n","Mesh \"aet3g2\" Dimension 2 Elemtype Triangle Nnode 3");
    fprintf(fp,"%s\n","Coordinates");
    fprintf(fp,"%s\n","End coordinates");
    fprintf(fp,"%s\n","Elements");
    for (i=1; i<=nelem[2]; ++i)
    {
        fprintf(fp,"%d",i);
        for (j=1; j<=nnode[2]; ++j)
        {
            fprintf(fp," %d",node[++n]);
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"%s\n","End elements");
    fprintf(fp,"%s\n","Mesh \"aeq9g3\" Dimension 2 Elemtype Quadrilateral Nnode 9");
    fprintf(fp,"%s\n","Coordinates");
    fprintf(fp,"%s\n","End coordinates");
    fprintf(fp,"%s\n","Elements");
    for (i=1; i<=nelem[3]; ++i)
    {
        fprintf(fp,"%d",i);
        for (j=1; j<=nnode[3]; ++j)
        {
            fprintf(fp," %d",node[++n]);
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"%s\n","End elements");
    fprintf(fp,"%s\n","Mesh \"aet6\" Dimension 2 Elemtype Triangle Nnode 6");
    fprintf(fp,"%s\n","Coordinates");
    fprintf(fp,"%s\n","End coordinates");
    fprintf(fp,"%s\n","Elements");
    for (i=1; i<=nelem[4]; ++i)
    {
        fprintf(fp,"%d",i);
        for (j=1; j<=nnode[4]; ++j)
        {
            fprintf(fp," %d",node[++n]);
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"%s\n","End elements");
    fclose(fp);
    if ((fp = fopen("./a1.gid/a1.post.res","w"))==NULL)
    {
        ;
        printf("cat:cannot open %s\n", "a1.post.res");
        return  1;
    }
    fprintf(fp,"%s\n","GID Post Results File 1.0");
    fclose(fp);
    return 0;
}
