function[out]=expand(U)
    [M,N]=size(U); %get input field array size
    K=256/M;
    % RM=rand([M/16,N/16])*2*pi;
    out = kron(U,ones(K));
end