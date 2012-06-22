%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : K-SVD algorithm for learning a dictionary
%-------------------------------------------------------------------------
% Inputs
%   X     : Data as columns
%   L     : Number of dictionary atoms
%   niter : Number of K-SVD iterations
%
% Ouputs
%   U     : Learned dictionary
%
function [U, info] = learn_dict(X, L, niter)
  if L < 2
    U = [];
    info = [];
    return;
  end

  % Initialize dictionary with random data patches
  idx = randsample(size(X, 2), L);
  U = X(:, idx);
  
  % Remove means and normalize
  U = U - repmat(mean(U), size(U, 1), 1);
  U = U ./ repmat(sqrt(sum(U.^2)), size(U, 1), 1);
  
  % Set first element to constant (the only one with non-zero mean)
  U(:, 1) = ones(size(U(:, 1)));
  U(:, 1) = U(:, 1) / norm(U(:, 1));
  
  % Mask to be used in sparse coding (use all data elements)
  M = ones(size(X));

  % K-SVD loop
  for i = 1:niter
    disp(['* Iteration: ' num2str(i)]);
    if i > 1
      U = remove_duplicates(U, Z, X);
    end
    % Sparsely decompose data according to current dictionary
    Z = omp(U, X, M, 0.05);
    print_info(U, Z, X);
    for j = randperm(L)
      % Don't update constant term
      if j == 1
        continue;
      end
      % Only consider nonzero elements
      Zm = Z;
      Zm(j, :) = 0;
      idx = Z(j, :) ~= 0;
      E = X - U*Zm;
      E = E(:, idx);
      % Find closest rank-1 approximation
      [u1, s1, z1] = svds(E, 1);
      % Remove mean from new element and renormalize
      u1 = u1 - mean(u1);
      u1 = u1 / norm(u1);
      % Update dictionary and sparse code
      U(:, j) = u1;
      Z(j, idx) = s1*z1;
    end
    print_info(U, Z, X);
    %figure; plot_dict(U); drawnow;
  end
  info.nnz = nnz(Z)/numel(Z);
end

function print_info(U, Z, X)
  display(['residual = ', num2str(norm(X - U*Z, 'fro'))]);
  display(['nnz = ', num2str(100*nnz(Z)/numel(Z)), '%']);
end

% Remove elements that are very similar and elements that are not used much
function U = remove_duplicates(U, Z, X)
  Unew = U;
  for i = 2:size(U, 2)
    Ui = U;
    Ui(:, i) = [];
    tmp = Ui'*U(:, i);
    if max(tmp) > 0.95 || nnz(Z(i, :)) < 5
      display(['removing ', num2str(i)]);
      Unew(:, i) = X(:, randsample(size(X, 2), 1));
      Unew(:, i) = U(:, i) - mean(U(:, i));
      Unew(:, i) = U(:, i) / norm(U(:, 1));
    end
  end
  U = Unew;
end