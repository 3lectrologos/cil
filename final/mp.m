%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : Modified matching pursuit for dealing with missing pixels
%-------------------------------------------------------------------------
% Inputs
%   U      : Dictionary with atoms as columns
%   X      : Masked data as columns
%   M      : Mask of unknown data
%   rc_min : Minimum allowed residual before stopping
%   sigma  : Target residual norm percentage
%
% Ouputs
%   Z      : Sparse code
%
function Z = mp(U, X, M, sigma, rc_min)
  l = size(U, 2);
  n = size(X, 2);
  Z = zeros(l, n);
  for i = 1:n
    x = X(M(:, i)==1, i);
    Ucur = U(M(:, i)==1, :);
    nf = sqrt(sum(Ucur.^2));
    Ucur = Ucur ./ repmat(nf, size(Ucur, 1), 1);
    residual = x;
    rc_max = Inf;
    while norm(residual) > sigma*norm(x) && rc_max > rc_min
      [rc_max, idx] = max(abs(Ucur'*residual));
      idx = idx(1);
      u = Ucur(:, idx);
      prod = u'*residual;
      Z(idx, i) = Z(idx, i) + prod;
      residual = (residual - prod*u);
    end
    Z(:, i) = Z(:, i) ./ nf';
  end
end