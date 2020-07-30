function dx = Lorenz(x,par)

[dim,nop] = size(x);
%if size(par,1)~=6
%par=par';end

dx = zeros(dim,nop);
dx(1,:) = par(1,:).*(x(2,:)-x(1,:));
dx(2,:) = par(2,:).*x(1,:) - x(2,:) -x(1,:).*x(3,:);
dx(3,:) = -par(3,:).*x(3,:) + x(1,:).*x(2,:);

%function y = lorenz (x, t)
%y = [10*(x(2) - x(1));
%x(1)*(28 - x(3)) - x(2);
%x(1)*x(2) - 8/3*x(3)];
%endfunction