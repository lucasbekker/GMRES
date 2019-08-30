function [ x, r, i ] = GMRES ( A, b, x0, tol, maxit )

r0 = A*x0 - b;
beta = norm(r0,2);
V(:,1) = r0/beta;
g = beta;
Q = cast(1,'like',b);

i = 0;
while ((abs(g(i + 1)) > tol) && (i < maxit))
    % Arnoldi %
    i = i + 1;
    w = A*V(:,i);
    for j = 1:i
        H(j,i) = w'*V(:,j);
        w = w - H(j,i)*V(:,j);
    end
    H((i + 1),i) = norm(w,2);
    V(:,(i + 1)) = w/H((i + 1),i);
    
    % Construct R and Givens rotation %
    H(1:i,i) = Q*H(1:i,i);
    rho = H(i,i);
    H(i,i) = sqrt(rho^2 + H((i + 1),i)^2);
    c = rho/H(i,i);
    s = H((i + 1),i)/H(i,i);
    H((i + 1),i) = 0;
    
    % Apply Givens rotation to Q %
    Q((i + 1),:) = -s*Q(i,:);
    Q(i,:) = c*Q(i,:);
    Q((i + 1),(i + 1)) = c;
    Q(i,(i + 1)) = s;
    
    % Apply Givens rotation to g %
    g(i + 1,1) = -s*g(i,1);
    g(i,1) = c*g(i,1);
end

y = H((1:i),:)\g(1:i);

x = x0 - V(:,(1:i))*y;
r = abs(g(i + 1));

end