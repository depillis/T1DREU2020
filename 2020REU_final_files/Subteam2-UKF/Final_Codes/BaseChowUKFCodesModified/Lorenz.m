% This function describes the Lorenz system

function dx = Lorenz(x,par)

[dim,nop] = size(x);
dx = zeros(dim,nop); %Initialize matrix for differential equations


%Parameters are scalar constants that are multiplied by the x values
dx(1,:) = par(1)*(x(2,:)-x(1,:)); %equation 48 from Chow
dx(2,:) = par(2)*x(1,:) - x(2,:) -x(1,:).*x(3,:); %equation 49 from Chow
dx(3,:) = -par(3)*x(3,:) + x(1,:).*x(2,:); %equation 50 from Chow

