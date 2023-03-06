function [abc, Vcol]=IsoTet4ShapeFunDer(xyztet)


x1=xyztet(1,1);x2=xyztet(1,2);x3=xyztet(1,3);x4=xyztet(1,4);
y1=xyztet(2,1);y2=xyztet(2,2);y3=xyztet(2,3);y4=xyztet(2,4);
z1=xyztet(3,1);z2=xyztet(3,2);z3=xyztet(3,3);z4=xyztet(3,4);


x12=x1-x2; x13=x1-x3; x14=x1-x4; x23=x2-x3; x24=x2-x4; x34=x3-x4;
x21=-x12; x31=-x13; x41=-x14; x32=-x23; x42=-x24; x43=-x34;
y12=y1-y2; y13=y1-y3; y14=y1-y4; y23=y2-y3; y24=y2-y4; y34=y3-y4;
y21=-y12; y31=-y13; y41=-y14; y32=-y23; y42=-y24; y43=-y34;
z12=z1-z2; z13=z1-z3; z14=z1-z4; z23=z2-z3; z24=z2-z4; z34=z3-z4;
z21=-z12; z31=-z13; z41=-z14; z32=-z23; z42=-z24; z43=-z34;

%compute the a b c
a1=y42*z32-y32*z42; b1=x32*z42-x42*z32; c1=x42*y32-x32*y42;
a2=y31*z43-y34*z13; b2=x43*z31-x13*z34; c2=x31*y43-x34*y13;
a3=y24*z14-y14*z24; b3=x14*z24-x24*z14; c3=x24*y14-x14*y24;
a4=y13*z21-y12*z31; b4=x21*z13-x31*z12; c4=x13*y21-x12*y31;
abc=[a1 b1 c1;
	a2 b2 c2;
	a3 b3 c3;
	a4 b4 c4];

% compute the volume
V01=x2*(y3*z4-y4*z3)+x3*(y4*z2-y2*z4)+x4*(y2*z3-y3*z2);
V02=x1*(y4*z3-y3*z4)+x3*(y1*z4-y4*z1)+x4*(y3*z1-y1*z3);
V03=x1*(y2*z4-y4*z2)+x2*(y4*z1-y1*z4)+x4*(y1*z2-y2*z1);
V04=x1*(y3*z2-y2*z3)+x2*(y1*z3-y3*z1)+x3*(y2*z1-y1*z2);
Vcol=(V01+V02+V03+V04)/6;

end