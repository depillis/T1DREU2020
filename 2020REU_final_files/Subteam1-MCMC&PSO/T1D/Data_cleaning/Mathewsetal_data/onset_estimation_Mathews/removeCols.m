function colstodelete =  removeCols(data)

colstodelete = zeros(1,1);
for i = 2:2:(length(data(1,:)))
    col = data(:, i);
    if ones(length(data(:,1)), 1) == isnan(col)
       colstodelete = [colstodelete i i-1];
    end
end
colstodelete(1) = [];
end