function breakpoints = findbreakpoints( x )
    n = length(x);
    diff1 = diff(x);
    breakpoints = zeros(n:1);
    breakpoints(1) = 1;
    count = 1;
    
    for i = 2:n-1
        q = diff1(i);
        w = diff1(i-1);
        if abs(q-w)>1e-4
            count = count+1;
            breakpoints(count) = i;
        end
    end
    count = count+1;
    breakpoints(count) = n;
    breakpoints = breakpoints(1:count);
end

