\ ------------------------------- equation --------------------------------------
\ For magnetic field in 2dxy(two dimensional cartesian coordinate system).
\ Compute magnetic field using least square method as
\                             H=1/fmu*curl(Az)
\ -------------------------------------------------------------------------------
DISP hx hy
COEF az
COOR x y
SHAP %1 %2
GAUS %3
MASS %1 vol
$CC double fmu,curlazx,curlazy;
MATE fmu 12.56637e-7

STIF
$CC vol = 1.0;
$CV curlazx=+{az/y}
$CV curlazy=-{az/x}
DIST=+[hx;hx]*0.0

LOAD=+[hx]*vol/fmu*curlazx+[hy]*vol/fmu*curlazy


END
