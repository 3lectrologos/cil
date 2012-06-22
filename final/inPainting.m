%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : Image inpainting pipeline combine Gaussian process
%               regression, intensity propagation, and sparse coding
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
  % Use dictionary?
  use_dict = true;
  % Use intensity propagation?
  use_ip = true;
  
  % Threshold parameters for choosing between methods
  marglik_thr = -0.8;
  sparsity_thr = 85;
  
  % GP function and hyperparameter setup
  gpp.meanfunc = @meanConst; gpp.hyp.mean = 0.5;
  gpp.covfunc = {@covMaterniso, 5}; gpp.hyp.cov = [1.2 -1.5];
  gpp.likfunc = @likGauss; gpp.hyp.lik = -3;
  
  % Initialize and split into blocks
  M = M ~= 0;
  d = 512 / k;
  dim = k * ones(1, d);
  CI = mat2cell(I, dim, dim);
  CM = mat2cell(M, dim, dim);
  CR = cell(d, d);
  
  % Load dictionary
  temp = load('dictionary.mat');
  U = temp.U;
  
  % Perform intensity propagation
  I_ec = ec(I, M); 
  CI_ec = mat2cell(I_ec, dim, dim);
  
  % Main loop combining methods for each patch
  for i = 1:size(CI, 1)
    for j = 1:size(CI, 2)
      Ib = CI{i, j};
      Mb = CM{i, j};
      if nnz(Mb) == numel(Mb)
        CR{i, j} = Ib;
        continue;
      end
      % GP training and test points
      [x1, x2] = meshgrid(1:k, 1:k);
      xtrain = [x2(Mb) x1(Mb)];
      ytrain = Ib(Mb);
      xtest = [x2(:) x1(:)];
      % Set GP mean prior to the mean of the current patch
      gpp.hyp.mean = mean(ytrain);
      % Compute GP log marginal likelihood
      lik = gp(gpp.hyp, @infExact, gpp.meanfunc, gpp.covfunc, gpp.likfunc, xtrain, ytrain);
      lik = lik / nnz(Mb);
      cflag = use_ip;
      % If the marginal likelihood is high enough use GP + IP
      if lik > marglik_thr
        % Else try sparse reconstruction
        [sr, nz] = sparse_rec(Ib, Mb, U, k);
        if use_dict && nz <= sparsity_thr
          cflag = false;
          CR{i, j} = sr;
        else
          % If the reconstruction is not sparse enough resort to GP + IP
          [CR{i, j}, v] = gp_rec(gpp, k, xtrain, ytrain, xtest);
        end
      else
        [CR{i, j}, v] = gp_rec(gpp, k, xtrain, ytrain, xtest);
      end
      % Combine GP and IP based on GP inferred variance
      if cflag
        V = min(20*reshape(v, k, k), 1);
        CR{i, j} = (1-V).*CR{i, j} + V.*CI_ec{i, j};
      end
      CR{i, j}(Mb) = Ib(Mb);
    end
  end
  R = cell2mat(CR);
  R(R > 1) = 1;
  R(R < 0) = 0;
end

% Reconstruct missing pixels using GP regression
function [R, v] = gp_rec(gpp, k, xtrain, ytrain, xtest)
  [m, v] = gp(gpp.hyp, @infExact, gpp.meanfunc, gpp.covfunc, gpp.likfunc, xtrain, ytrain, xtest);
  R = reshape(m, k, k);
end

% Reconstruct missing pixels using sparse coding
function [R, nz] = sparse_rec(Ib, Mb, U, k)
  Xr = reshape(Ib, [], 1);
  Mr = reshape(Mb, [], 1);
  Z = mp(U, Xr, Mr, 0.03, 0.01);
  R = reshape(U*Z, k, k);
  nz = nnz(Z);
end