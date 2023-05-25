function [out] = reducemat(U,a)
    [M,N]=size(U); %get input field array size
    K=M/a;
    out = ones(a,a);
    for ii = 1:a
        for jj = 1:a
            out(ii,jj) = out(ii,jj) * U(ii*16-1,jj*16-1)
        end
    end
end