function[out]=chk_J(Target,T_Size);
[M,N]=size(Target);
out=sum(sum(Target(M/2-T_Size/2+1 : M/2+T_Size/2, N/2-T_Size/2+1 : N/2+T_Size/2)));
end