function A = OffsetShift(A, mask0)

    a0 = A.*mask0;
    a0 = sum(a0(:)) ./ sum(mask0(:));
%     a0(a0==0) = [];

    A = A -mean(a0);
    A(A<0) = 0;

end