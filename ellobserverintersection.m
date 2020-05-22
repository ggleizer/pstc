function [z,Z] = ellobserverintersection(x,y,X,Y,C)
%ELLOBSERVERINTERSECTION Computes the trace-optimal ellipsoid for given
%predicted state in ell(x,X) and input with noise in ell(y,Y), for a given
%C that satisfies y = Cx.
%
%   [z,Z] = ELLOBSERVERINTERSECTION(x,y,X,Y,C)
%
%   Inputs:
%       x in R^n
%       y in R^m
%       X = X', X > 0, X in R^(n x n)
%       Y = Y', Y > 0, Y in R^(m x m)
%       C in R^(m x n)
%
%   Outputs:
%       z in R^n
%       Z = Z', Z >= 0, Z in R^(n x n).
%
%   If x_k in ell(x,X) and y_k in ell(y,Y), then x_k in ell(z,Z). This is
%   the trace-optimal ellipsoid that contains the intersection of the given
%   ellipsoids, among a family of parametrized ellipsoids. The optimization
%   is done via fminbnd on the sole parameter that needs to be optimized.
% 
%   Author: Gabriel de A. Gleizer, Aug 2018 (g.gleizer@tudelft.nl)

e = y - C*x;
Xi = inv(X);
Yi = inv(Y);

f = @(l) lambdafusiontrace(l,X,Y,Xi,Yi,C,e);

l = fminbnd(f,sqrt(eps),1);  % 0 < l <= 1

Z = l*Xi + (1-l)*C'*Yi*C;
a = 1 - l*(1-l)*e'*((l*Y + (1-l)*C*X*C')\e);
Zi = inv(Z);
Z = a*Zi;

z = Zi*(l*Xi*x + (1-l)*C'*Yi*y);

end

function tr = lambdafusiontrace(l,X,Y,Xi,Yi,C,e)
    Z = l*Xi + (1-l)*C'*Yi*C;
    a = 1 - l*(1-l)*e'*((l*Y + (1-l)*C*X*C')\e);
    tr = a*trace(inv(Z));
end