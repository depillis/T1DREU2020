function y = my_linspace(d1, d2, n)
% MY_LINSPACE: Generate n evenly spaced numeric array between array d1, d2
if nargin == 2
    n = 100; % Default number of nodes
end
n = double(round(n));
if all(size(d1) == size(d2))
    siz = size(d1);
    y = zeros([siz n]);
    n1 = floor(n)-1; n2 = zeros(1,1,n1+1);n2(:)= (0:n1);
    for i=1:numel(d1)
        [k,l]=ind2sub(siz,i);
        c = (d2(i) - d1(i))*(n1-1); % opposite signs may cause overflow
        if isinf(c)
            y(k,l,:) = d1(i) + (d2(i)/n1).*n2 - (d1(i)/n1).*n2; 
        else
            y(k,l,:) = d1(i) + n2.*(d2(i) - d1(i))/n1;
        end
        y(k,l,1) = d1(i);
        y(k,l,end) = d2(i);
    end
else
    error('Vectors must be same length')
end