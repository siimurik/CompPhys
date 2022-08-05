function y = f(x)
%    y = x .* sin (1./x) .* sqrt (abs (1 - x));
    y = x .* cos(10.0 .* x.^2) ./ (x.^2 + 1.0);
endfunction