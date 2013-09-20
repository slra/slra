% convert misfit into the System ID toolbox compare criterion

function c = m2c(M, u, y, eiv, cmp)

if cmp
    T = size(y, 1);
    if eiv
        nw = norm([u y] - ones(T, 1) * mean([u y]), 'fro');
    else
        nw = norm(y - ones(T, 1) * mean(y), 'fro');
    end
    c = 100 * max(0, 1 - sqrt(M) / nw); 
else
    if eiv
        nw = norm([u y], 'fro');
    else
        nw = norm(y, 'fro');
    end
    c = 100 * sqrt(M) / nw;
end