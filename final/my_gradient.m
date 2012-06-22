%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : Custom image gradient
%-------------------------------------------------------------------------
% Inputs
%   I    : Masked image
%   mask : Mask of missing pixels
%
% Ouputs
%   Ix   : Gradient in the x-direction
%   Iy   : Gradient in the y-direction
%
function [Ix, Iy] = my_gradient(I, mask)
  I(~mask) = nan;

  Ix1 = nan(size(I));
  Ix2 = nan(size(I));
  Iy1 = nan(size(I));
  Iy2 = nan(size(I));

  Ix1(:, 1:end - 1) = I(:, 2:end) - I(:, 1:end - 1);
  Ix2(:, 2:end) = I(:, 2:end) - I(:, 1:end - 1);
  Iy1(1:end - 1, :) = I(2:end, :) - I(1:end - 1, :);
  Iy2(2:end, :) = I(2:end, :) - I(1:end - 1, :);

  Mx1 = ~isnan(Ix1);
  Mx2 = ~isnan(Ix2);
  My1 = ~isnan(Iy1);
  My2 = ~isnan(Iy2);

  Ix1(isnan(Ix1)) = 0;
  Ix2(isnan(Ix2)) = 0;
  Iy1(isnan(Iy1)) = 0;
  Iy2(isnan(Iy2)) = 0;

  Ix = (Ix1 + Ix2) ./ (Mx1 + Mx2);
  Iy = (Iy1 + Iy2) ./ (My1 + My2);

  Ix(isnan(Ix) & mask) = 0;
  Iy(isnan(Iy) & mask) = 0;
end