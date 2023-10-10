function Quiver = QuiverShift_v2(h,Quiver)
% Quiver = Quiver_RI2;    
%%
    mask = ones(size(Quiver.U), 'single');
    mask(3:end-2, 3:end-2, 3:end) = 0;
    
    U0 = Quiver.U .* mask;
    V0 = Quiver.V .* mask;
    W0 = Quiver.W .* mask;

    U0 = U0(:);
    V0 = V0(:);
    W0 = W0(:);
    U0(~mask(:)) = [];
    V0(~mask(:)) = [];
    W0(~mask(:)) = [];

    I = Quiver.I; I = I - mean(I(:));
    J = Quiver.J; J = J - mean(J(:));
    K = Quiver.K; K = K - mean(K(:));

    I = I*10^-6;
    J = J*10^-6;
    K = K*10^-6;

    I0 = I .* mask;
    J0 = J .* mask;
    K0 = K .* mask;
    I0(~mask(:)) = [];
    J0(~mask(:)) = [];
    K0(~mask(:)) = [];

% prevent exceeding indices
        epsilon = 0;
        
% 
%         A0 = [I(:).^2 J(:).^2 I(:) J(:) ones(length(I(:)),1)];
%         A = [I0(:).^2 J0(:).^2 I0(:) J0(:) ones(length(I0(:)),1)];
% 
%         y1=(A'*A+epsilon)\(A'*U0(:));
%         Quiver.U(:) = Quiver.U(:) - A0*y1;
% 
%         y2=(A'*A+epsilon)\(A'*V0(:));
%         Quiver.V(:) = Quiver.V(:) - A0*y2;

%         A = [K0(:).^2 K0(:) ones(length(I0(:)),1)];
%         y3=(A'*A+epsilon)\(A'*W0(:));
%         A0 = [K(:).^2 K(:) ones(length(I(:)),1)];
%         Quiver.W(:) = Quiver.W(:) - A0*y3;


        A = [I0(:).^2 J0(:).^2 K0(:).^2 I0(:).*J0(:) I0(:).*K0(:) I0(:) J0(:) K0(:) ones(length(I0(:)),1)];
        y1=(A'*A+epsilon)\(A'*U0(:));
        A = [I(:).^2 J(:).^2 K(:).^2 I(:).*J(:) I(:).*K(:) I(:) J(:) K(:) ones(length(I(:)),1)];
        Quiver.U(:) = Quiver.U(:) - A*y1;

        A = [I0(:).^2 J0(:).^2 K0(:).^2 J0(:).*I0(:) J0(:).*K0(:) I0(:) J0(:) K0(:) ones(length(I0(:)),1)];
        y2=(A'*A+epsilon)\(A'*V0(:));
        A = [I(:).^2 J(:).^2 K(:).^2 J(:).*I(:) J(:).*K(:) I(:) J(:) K(:) ones(length(I(:)),1)];
        Quiver.V(:) = Quiver.V(:) - A*y2;

        A = [I0(:).^2 J0(:).^2 K0(:).^2 K0(:).*I0(:) J0(:).*K0(:) I0(:) J0(:) K0(:) ones(length(I0(:)),1)];
        y3=(A'*A+epsilon)\(A'*W0(:));
        A = [I(:).^2 J(:).^2 K(:).^2 K(:).*I(:) J(:).*K(:) I(:) J(:) K(:) ones(length(I(:)),1)];
        Quiver.W(:) = Quiver.W(:) - A*y3;

end