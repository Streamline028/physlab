function[out]=radCal_zero(U)
[M,N]=size(U); %get input field array size
for i=M/2:M
    if (U(i,N/2)<=0)
        out = i - M/2;
        break
    end
end
end