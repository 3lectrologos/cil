% Measure approximation error and compression ratio for several images.
% NOTE: Images must be have .png ending.
function EvaluateInpainting(varargin)
  NUM_MASKS = 2;
  ROOT_DIR = fileparts(fileparts(mfilename('fullpath')));
  DATA_DIR = fullfile(ROOT_DIR, 'data');
  file_list = dir(DATA_DIR);
  if nargin > 0
    show = varargin{1};
  else
    show = false;
  end

  tot = [];
  for j = 1:NUM_MASKS
    mask_name = sprintf('mask_%d.png', j);
    errors = [];
    k = 1;
    display(['Using ', mask_name]);
    display('----------------');
    for i = 3:length(file_list)
      file_name = file_list(i).name;
      % Only keep the images in the loop
      if (length(file_name) < 5)
        continue;
      elseif (max(file_name(end-6:end) ~= '512.png'))
        continue;
      end
      % Read image, convert to double precision and map to [0,1] interval
      I = imread(fullfile(DATA_DIR, file_name));
      I = double(I) / 255;
      % Read the respective binary mask
      % EVALUATION IS DONE WITH A FIXED MASK
      mask = imread(fullfile(DATA_DIR, mask_name));
      I_mask = I;
      I_mask(~mask) = 0;
      % Call the main inPainting function
      I_rec = inPainting(I_mask, mask);
      % Measure approximation error
      error = mean(mean(mean((I - I_rec).^2)));
      errors = [errors error];
      fprintf(1, '%s: %f\n', file_name, error);
      k = k+1;
      % Show result
      if show
        imsize = [1200 500];
        screensize = get(0, 'ScreenSize');
        xpos = ceil((screensize(3)-imsize(1))/2);
        ypos = ceil((screensize(4)-imsize(2))/2);
        figure('Position', [xpos ypos imsize]);
        subplot(1, 3, 1); imshow(I); title('Original');
        subplot(1, 3, 2); imshow(I_mask); title('Masked');
        I_rec(I_rec > 1) = 1;
        I_rec(I_rec < 0) = 0;
        subplot(1, 3, 3); imshow(I_rec); title('Reconstructed');
        spaceplots([0.01 0.01 0.01 0.01], [0.01 0.01]);
        drawnow;
      end
    end
    tot = [tot mean(errors)];
    disp(['Average quadratic error: ' num2str(mean(errors))]);
  end
  display('----------------------------');
  disp(['Total average quadratic error: ' num2str(mean(tot))]);
end