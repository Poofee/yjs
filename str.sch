\ ��ʽ������֪�������������ĳ�ֱ��ʽ��ֵ�ֲ�����С���˷�str.sch
\ �ռ���ɢ�������ʽΪ��
\           [S][U]=[F]
\ ---------------------------------------------------
DEFI
STIF s
MASS m
LOAD f
TYPE e
MDTY l

EQUATION
\............ ���Է����������(���о���) .............../
L,MASS = [m]
\................. ���Է������Ҷ��� ..................../
FORC = [f]

SOLUTION w
VECT w
$CC if (init == 0) init = 1;
$cc // �洢����w
WRITE(s,unod) w

@NBDE
 ntype = nbdetype;


END

