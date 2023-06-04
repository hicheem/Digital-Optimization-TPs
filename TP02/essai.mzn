var int: x;
var int: y;
%var int: f=150*x+450*y;
var int: f; constraint f=150*x+450*y;
constraint x>=0 /\ y>=0 /\ x<=120 /\ y<=70 /\ x+y<=140 /\ x+2*y <=180;
%var int: f; constraint f=160*x+300*y;
%constraint x>=0 /\ y>=0 /\ x<=120 /\ y<=70 /\ x+y <=140 /\ 2*x+3*y<=181;

solve maximize f;

output ["x=\(x);\ny=\(y);" ];
%output ["x=\(x);"]++["\ny=\(y);" ];


