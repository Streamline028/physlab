function[out]=radCal(U)
[M,N]=size(U); %get input field array size
for i=M/2:M
    if (U(i,N/2)<(max(max(U))*0.135))
        out = i - M/2;
        break
    else
        out = -90;
    end
end
end