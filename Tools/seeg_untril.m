function M = seeg_untril(V)

c = length(V);
n = (2*c + 1/4)^(1/2) + 1/2;

id = [];
for ci = 1:n-1       % Go through columns
for ri = ci+1:n      % Each column select subset of rows
    id = [id; ri, ci];
end
end

M = zeros(n);
for v = 1:c
    M(id(v,1), id(v,2)) = V(v); 
    M(id(v,2), id(v,1)) = V(v);
end
