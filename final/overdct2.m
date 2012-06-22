%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : Create a two-dimensional overcomplete DCT dictionary
%-------------------------------------------------------------------------
% Inputs
%   n : Square root of dictionary atom dimension
%   L : Number of dictionary atoms to create
%
% Ouputs
%   U : Return dictionary with atoms as columns
%
function U = overdct2(n, L)
  r = ceil(sqrt(L));
  U = zeros(n, r);
  for k = 0:1:r-1
    v = cos((0:1:n-1)'*k*pi/r);
    if k > 0
      v = v - mean(v);
    end;
    U(:, k+1) = v / norm(v);
  end;
  U = kron(U, U);
end