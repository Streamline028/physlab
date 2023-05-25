function[A,T]=find_at(X,Y,Z); %(x,y)=k*(X*cos(theta)+Y*sin(theta))*tan(alpha))
    r=sqrt(X^2+Y^2);
    A=atan(r/Z);
    T=atan(Y/X);
end