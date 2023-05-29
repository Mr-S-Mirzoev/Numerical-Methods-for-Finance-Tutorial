function error = heateq_fdm(a,b,T,N,M,theta)
h = (b-a)/(N+1); k = T\M; e = ones(N,1);
I = speye(N); G = h^(-2)*spdiags([-e 2*e -e],-1:1,N,N);
B = I+k*theta*G; C = I-k*(1-theta)*G
x = [a + h : h : b - h]'; u0 = x.*sin(pi * x)
f = -(1-pi^2)*x.*sin(pi*x)-2*pi*cos(pi*x);
u = zeros(N,T/k+1); u(:,1) = u0;
for j = 1:M
F = f*(theta*exp(-j*k)+(1-theta)*exp(-(j-1)*k));
u(:,j+1) = B\(C*u(:,j)+k*F);
err(j) = norm(u(:,j+1)-exp(-k*j)*u0);
end
error = sqrt(h)*max(err)