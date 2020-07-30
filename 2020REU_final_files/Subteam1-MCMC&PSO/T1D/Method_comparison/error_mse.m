% Mean squared error
function mse = error_mse(observed, predicted)% Sum-of-squares function
n = length(observed); 
ss = sum((observed-predicted).^2);

mse = ss/length(observed);



end