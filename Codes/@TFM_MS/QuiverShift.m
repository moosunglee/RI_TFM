function [Quiver, shift_vector] = QuiverShift(h,Quiver)
    mask = ones(size(Quiver.U), 'single');
    mask(3:end-2, 3:end-2, 3:end-2) = 0;
    
    U = Quiver.U .* mask;
    V = Quiver.V .* mask;
    W = Quiver.W .* mask;

    U = U(:);
    V = V(:);
    W = W(:);
    U(~mask(:)) = [];
    V(~mask(:)) = [];
    W(~mask(:)) = [];
    
    U(isoutlier(U)) = [];
    V(isoutlier(V)) = [];
    W(isoutlier(W)) = [];

    shift_vector = [mean(U) mean(V) mean(W)];

    Quiver.U = Quiver.U - shift_vector(1);
    Quiver.V = Quiver.V - shift_vector(2);
    Quiver.W = Quiver.W - shift_vector(3);
    
end