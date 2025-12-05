Variables x2, x1, objvar ;
Equations obj, c1, c2 ;

obj.. objvar =E= 0 ;
c1.. -1 + x2 =E= 0 ;
c2.. x1 - x2 =G= 0 ;

x2.lo = 0 ; 
x2.up = 1000 ;
x1.lo = 0 ; 
x1.up = 2000 ;
objvar.lo = -100 ; 
objvar.up = 100 ;
