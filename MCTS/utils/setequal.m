function e = setequal(a, b)
    e = isempty(setdiff(a,b)) && isempty(setdiff(b,a));
end