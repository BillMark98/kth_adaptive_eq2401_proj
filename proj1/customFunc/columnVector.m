function v = columnVector(v)

% given a vector v, convert it to standard vector
% form so of shape n * 1

[m,n] = size(v);
if (n ~= 1)
    if (m ~= 1)
        error("not a vector given")
    end
    v = v';
end