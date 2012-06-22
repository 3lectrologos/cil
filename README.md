A Combined Scheme for Image Inpainting
============================================
Final project of the _Computational Intelligence Lab_
(ETH Zurich, Spring Semester 2012).

Authors
-------
Alkis Gkotovos (<alkisg@student.ethz.ch>),
Bo Li (<libo@student.ethz.ch>) and
Alexey Sizov (<asizov@student.ethz.ch>)

Description
-----------
A method for image inpainting combining Gaussian Process regression,
intensity propagation, and sparse coding via an overcomplete dictionary.

Files
-----
The directory contains the files used by our combined method, as well as
the files of the baseline methods used for comparison. Note that the
main function for all methods is called `inPainting.m` and all of them
adhere to the same interface.

### Combined method
The files of our combined method reside under the `final/` subdirectory.
Included are also the files `covMaterniso.m`, `meanConst.m`, `likGauss.m`,
`infExact.m`, `solve_chol`, and `sq_dist.m` from the GPML toolbox
(http://www.gaussianprocess.org/gpml/code/matlab/doc/) that are used for
doing Gaussian process inference.

### GP baseline
The GP baseline can be run using the code of the combined method, after
setting the flags `use_dict` and `use_ip` in file `inPainting.m` to `false`.

### Dictionary baseline
The files of the dictionary baseline method can be found under the
`baselines/dict/` subdirectory.

### SVD baseline
The files of the SVD baseline method can be found under the `baselines/svd/`
subdirectory.

Usage
-----
For evaluating any of the methods on our image test set, the following steps
should be taken:

* Execute in Matlab the `startup.m` script found at the top directory.
* Change the working directory to one of the methods, i.e. to
  `final/` for running the combined method or the GP baseline,
  to `baselines/dictionary/` for running the dictionary baseline,
  or to `baselines/svd/` for running the svd baseline.
* Execute in Matlab the `EvaluateInpainting.m` script, which should have
  already been added to the path after the first step, to run the
  corresponding method on all test images available. Providing a
  `true` argument to `EvaluateInpainting` will additionally print the
  reconstruction results as they are computed.

As an example, to evaluate our combined method and display the resulting
reconstructed images, use the following commands in Matlab:
```
> startup
> cd final
> EvaluateInpainting(true)
```