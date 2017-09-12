*npoin  *nelem
-1000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*loop nodes
*format "%6i %12.6e %12.6e %12.6e"
  *NodesNum *NodesCoord(1) *NodesCoord(2)
*end
-2001 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*Set Cond volume-a1a *nodes
*Add Cond surface-a1a *nodes
*Add Cond line-a1a *nodes
*Add Cond point-a1a *nodes
*loop nodes *OnlyInCond
    *NodesNum *cond(1)
*end
-2001 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*Set Cond volume-a1a *nodes
*Add Cond surface-a1a *nodes
*Add Cond line-a1a *nodes
*Add Cond point-a1a *nodes
*loop nodes *OnlyInCond
    *NodesNum *cond(2)
*end
-2002 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*Set Cond volume-a1b *nodes
*Add Cond surface-a1b *nodes
*Add Cond line-a1b *nodes
*Add Cond point-a1b *nodes
*loop nodes *OnlyInCond
    *NodesNum *cond(1) *cond(3)
*end
-2002 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*Set Cond volume-a1b *nodes
*Add Cond surface-a1b *nodes
*Add Cond line-a1b *nodes
*Add Cond point-a1b *nodes
*loop nodes *OnlyInCond
    *NodesNum *cond(2) *cond(4)
*end
-4000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
-4 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*set cond Surface-aeq4g2 *elems
*loop elems *OnlyIncond
*format "%10i %10i %10i %10i %10i %10i "
*ElemsNum *elemsConec *cond(1)
*end
-3 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*set cond Surface-aet3g2 *elems
*loop elems *OnlyIncond
*format "%10i %10i %10i %10i %10i "
*ElemsNum *elemsConec *cond(1)
*end
-9 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*set cond Surface-aeq9g3 *elems
*loop elems *OnlyIncond
*format "%10i %10i %10i %10i %10i %10i %10i %10i %10i %10i %10i "
*ElemsNum *elemsConec *cond(1)
*end
-6 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*set cond Surface-aet6 *elems
*loop elems *OnlyIncond
*format "%10i %10i %10i %10i %10i %10i %10i %10i "
*ElemsNum *elemsConec *cond(1)
*end
-4 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*set cond Surface-beq4g2 *elems
*loop elems *OnlyIncond
*format "%10i %10i %10i %10i %10i %10i "
*ElemsNum *elemsConec *cond(1)
*end
-3 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*set cond Surface-bet3g2 *elems
*loop elems *OnlyIncond
*format "%10i %10i %10i %10i %10i "
*ElemsNum *elemsConec *cond(1)
*end
-9 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*set cond Surface-beq9g3 *elems
*loop elems *OnlyIncond
*format "%10i %10i %10i %10i %10i %10i %10i %10i %10i %10i %10i "
*ElemsNum *elemsConec *cond(1)
*end
-6 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
*set cond Surface-bet6 *elems
*loop elems *OnlyIncond
*format "%10i %10i %10i %10i %10i %10i %10i %10i "
*ElemsNum *elemsConec *cond(1)
*end
-5000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
