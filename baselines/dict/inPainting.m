%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : Image inpainting baseline using sparse coding via
%               overcomplete DCT dictionary and matching pursuit
%-------------------------------------------------------------------------
% Inputs
%   I     : Masked image
%   mask  : Mask of missing pixels
%
% Ouputs
%   I_rec : Reconstructed image
%
function I_rec = inPainting(I, mask)
  % Patch size
  k = 16;
  
  % Matching pursuit parameters
  rc_min = 0.01;
  sigma = 0.01;

  % Convert mask to logical
  mask = mask ~= 0;
  
  % Split image and mask to patch columns
  X = im2col(I, [k k], 'distinct');
  M = im2col(mask, [k k], 'distinct');

  % Create overcomplete DCT dictionary
  U = overdct2(k, 900);

  % Sparse coding via matching pursuite
  Z = mp(U, X, M, sigma, rc_min);

  % Reconstruct
  rec = col2im(U*Z, [k k], size(I), 'distinct');
  I_rec = mask.*I + (~mask).*rec;
  I_rec(I_rec > 1) = 1;
  I_rec(I_rec < 0) = 0;
end