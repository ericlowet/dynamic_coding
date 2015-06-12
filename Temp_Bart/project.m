function u1=project(v1,v2)

v1=v1(:);
v2=v2(:);

numerator=v1.'*v2;
denominator=v2.'*v2;

u1=numerator/denominator*v2;