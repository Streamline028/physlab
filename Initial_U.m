function[out]=Initial_U(X); %U=del_u matrix J=delJ=J(u+del_u)-J(u)
[M,N]=size(X);
out=rand([M,N])*2*pi;
end