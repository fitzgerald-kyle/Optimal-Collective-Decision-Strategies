function val = simpson(x, f)
    % assumes evenly-spaced subintervals
    hx = x(2)-x(1);
    val = f(1) + f(end);
    for i=2:length(f)-1
        if mod(i,2)
            val = val + 2*f(i);
        else
            val = val + 4*f(i);
        end
    end
    val = hx*val/3;
end