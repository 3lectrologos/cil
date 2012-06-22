%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : Orthogonal matching pursuit
%-------------------------------------------------------------------------
% Inputs
%   U        : Dictionary with atoms as columns
%   X        : Masked data as columns
%   M        : Mask of unknown data
%   sigma    : Target residual norm percentage
%   l (opt.) : Maximum number of basis elements (default: size of dict)
%
% Ouputs
%   Z        : Sparse code
%
function Z = omp(U, X, M, sigma, varargin)
  l = size(U, 2);
  n = size(X, 2);
  Z = zeros(l, n);
  if nargin > 4
    maxiter = ceil(l*varargin{1});
  else
    maxiter = l;
  end
  % Loop over all observations in the columns of X
  for i = 1:n
    x = X(M(:, i)==1, i);
    Ucur = U(M(:, i)==1, :);
    residual = x;
    L = [];
    idxs = [];
    iter = 1;
    while norm(residual) > sigma*norm(x) && iter < maxiter
      [~, idx] = max(abs(Ucur'*residual));
      idx = idx(1);
      idxs = [idxs idx];
      L = [L Ucur(:, idx)];
      [Q, ~] = qr(L, 0);
      residual = (x - Q*Q'*x);
      iter = iter + 1;
    end
    % Add the calculated coefficient vector z to the overall matrix Z
    Z(idxs, i) = pinv(L)*x;
  end
end