%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : Image inpainting baseline using SVD
%-------------------------------------------------------------------------
% Inputs
%   I : Masked image
%   M : Mask of missing pixels
%
% Ouputs
%   R : Reconstructed image
%
function R = inPainting(I, M)
  % Patch size
  k = 16;
  
  % Initialize and split into blocks
  M = M ~= 0;
  d = 512 / k;
  dim = k * ones(1, d);
  CI = mat2cell(I, dim, dim);
  CM = mat2cell(M, dim, dim);
  CR = cell(d, d);
  
  % Loop through patches
  for i = 1:size(CI, 1)
    for j = 1:size(CI, 2)
      Ib = CI{i, j};
      Mb = CM{i, j};
      if nnz(Mb) == numel(Mb)
        CR{i, j} = Ib;
        continue;
      end
      F = impute(Ib, Mb);
      [U, D, V] = svd(F);
      % Keep a percentage of the total variance
      c = find(cumsum(diag(D))/sum(diag(D)) > 0.75);
      k = c(1);
      CR{i, j} = U(:, 1:k)*D(1:k, 1:k)*V(:, 1:k)';
      CR{i, j}(Mb) = Ib(Mb);
    end
  end
  R = cell2mat(CR);
end

% Median-based imputation of missing values
function R = impute(I, M)
  R = I;
  NM = ~M;
  pred = zeros(1, size(R, 1));
  for i = 1:size(R, 1)
    buf = R(i, :);
    pred(i) = median(buf(M(i, :)));
    if(isnan(pred(i)))
      pred(i) = 0.5;
    end
  end
  predc = zeros(1, size(R, 2));
  for i = 1:size(R, 2)
    buf = R(:, i);
    predc(i) = median(buf(M(:, i)));
    if(isnan(predc(i)))
      predc(i) = 0.5;
    end            
  end
  X2 = R;
  for i = 1:size(R, 1)
    buf = R(i, :);
    buf(NM(i, :)) = pred(i);
    X2(i, :) = buf;
  end
  R = X2;
  for i = 1:size(R, 2)
    buf = R(:, i);
    buf(NM(:, i)) = 0.5*(buf(NM(:, i)) + predc(i));
    R(:, i) = buf;
  end
end