%-------------------------------------------------------------------------
% ETH Zurich, Spring Semester 2012
% Computation Intelligence Lab: Final project
%
% Authors     : Alkis Gkotovos <alkisg@student.ethz.ch>,
%               Bo Li <libo@student.ethz.ch> and
%               Alexey Sizov <asizov@student.ethz.ch>
% Description : Reconstruct missing image pixels using gradient-based
%               intensity propagation
%-------------------------------------------------------------------------
% Inputs
%   I     : Masked image
%   mask  : Mask of missing pixels
%
% Ouputs
%   I_rec : Reconstructed image
%
function I_rec = ec(I, mask)
  rec = fill_patch(I, mask);
  rec = medfilt2(rec);
  rec = medfilt2(rec);
  rec = imfilter(rec, fspecial('gaussian'));
  I = mask.*I + (~mask).*rec;
  I(I > 1) = 1;
  I(I < 0) = 0;
  I_rec = I;
end

function x = fill_patch(x, m)
  [dx, dy] = my_gradient(x, m);
  while ~all(m(:))
    x_est1 = nan(size(x));
    dx_est1 = nan(size(x));
    dy_est1 = nan(size(x));
    steep1 = inf(size(x));
    x_est1(:, 1:end - 1) = x(:, 2:end) - dx(:, 2:end);
    steep1(:, 1:end - 1) = abs(dx(:, 2:end));
    m1 = ~isnan(x_est1);
    x_est1(isnan(x_est1)) = 0;
    dx_est1(:, 1:end - 1) = dx(:, 2:end);
    dy_est1(:, 1:end - 1) = dy(:, 2:end);
    dx_est1(isnan(dx_est1)) = 0;
    dy_est1(isnan(dy_est1)) = 0;

    x_est2 = nan(size(x));
    dx_est2 = nan(size(x));
    dy_est2 = nan(size(x));
    steep2 = inf(size(x));
    x_est2(:, 2:end) = x(:, 1:end - 1) + dx(:, 1:end - 1);
    steep2(:, 2:end) = abs(dx(:, 1:end - 1));
    m2 = ~isnan(x_est2);
    x_est2(isnan(x_est2)) = 0;
    dx_est2(:, 2:end) = dx(:, 1:end - 1);
    dy_est2(:, 2:end) = dy(:, 1:end - 1);
    dx_est2(isnan(dx_est2)) = 0;
    dy_est2(isnan(dy_est2)) = 0;

    x_est3 = nan(size(x));
    dx_est3 = nan(size(x));
    dy_est3 = nan(size(x));
    steep3 = inf(size(x));
    x_est3(1:end - 1, :) = x(2:end, :) - dy(2:end, :);
    steep3(1:end - 1, :) = abs(dy(2:end, :));
    m3 = ~isnan(x_est3);
    x_est3(isnan(x_est3)) = 0;
    dx_est3(1:end - 1, :) = dx(2:end, :);
    dy_est3(1:end - 1, :) = dy(2:end, :);
    dx_est3(isnan(dx_est3)) = 0;
    dy_est3(isnan(dy_est3)) = 0;

    x_est4 = nan(size(x));
    dx_est4 = nan(size(x));
    dy_est4 = nan(size(x));
    steep4 = inf(size(x));
    x_est4(2:end, :) = x(1:end - 1, :) + dy(1:end-1, :);
    steep4(2:end, :) = abs(dy(1:end - 1, :));
    m4 = ~isnan(x_est4);
    x_est4(isnan(x_est4)) = 0;
    dx_est4(2:end, :) = dx(1:end - 1, :);
    dy_est4(2:end, :) = dy(1:end - 1, :);
    dx_est4(isnan(dx_est4)) = 0;
    dy_est4(isnan(dy_est4)) = 0;

    x_est = x;
    dx_est = dx;
    dy_est = dy;
    steep = inf(size(x));
    x_est(steep > steep1) = x_est1(steep > steep1);
    dx_est(steep > steep1) = dx_est1(steep > steep1);
    dy_est(steep > steep1) = dy_est1(steep > steep1);
    steep(steep > steep1) = steep1(steep > steep1);
    x_est(steep > steep2) = x_est2(steep > steep2);
    dx_est(steep > steep2) = dx_est2(steep > steep2);
    dy_est(steep > steep2) = dy_est2(steep > steep2);
    steep(steep > steep2) = steep2(steep > steep2);
    x_est(steep > steep3) = x_est3(steep > steep3);
    dx_est(steep > steep3) = dx_est3(steep > steep3);
    dy_est(steep > steep3) = dy_est3(steep > steep3);
    steep(steep > steep3) = steep3(steep > steep3);
    x_est(steep > steep4) = x_est4(steep > steep4);
    dx_est(steep > steep4) = dx_est4(steep > steep4);
    dy_est(steep > steep4) = dy_est4(steep > steep4);

    x(~m) = x_est(~m);
    dx(~m) = dx_est(~m);
    dy(~m) = dy_est(~m);

    m = m | (m1 | m2 | m3 | m4);
  end
end