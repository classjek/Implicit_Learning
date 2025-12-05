Variables x2, x1, objvar ;
Equations obj, c1, c2, c3, c4, c5, c6, c7, c8, c9 ;

obj.. objvar =E= 0 ;
c1..  - x1 + x1^2 =E= 0 ;
c2.. x1 - x1  =E= 0 ;
c3..  - x2 + x2^2 =E= 0 ;
c4.. x2 =E= 0 ;
c5.. -1.000000 + x2 =E= 0 ;
c6.. 1.000000 + x1 =G= 0 ;
c7.. 2.000000 - x1 =G= 0 ;
c8.. 1.000000 + x2 =G= 0 ;
c9.. 2.000000 - x2 =G= 0 ;

x2.lo = 0 ;
x2.up = 1000 ;
x1.lo = 0 ;
x1.up = 1000 ;
objvar.lo = -100;
objvar.up = 100;
