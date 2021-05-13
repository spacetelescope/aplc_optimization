The input file parameters (telescope aperture and Lyot stop specifications) are particular to the instrument being surveyed,
otherwise details on each user-definable survey parameter are as follows:

**Lyot Stop parameters**
 - ``alignment_tolerance``: Lyot stop alignment tolerance, in pixels.
 - ``num_lyot_stops``: number of Lyot stop translations.

**Focal Plane Mask parameters**
 - ``radius``: radius of the FPM of the APLC, in lambda_0/D.
 - ``num_pix``: number of pixels in the FPM.
 - ``grayscale```: grayscale FPM, otherwise black & white.
 - ``field_stop_radius``: radius of the FPM stop.

**Image parameters**
 - ``contrast``: contrast goal in the dark zone of the coronagraphic image.
 - ``iwa``: inner edge of dark zone region, in lambda_0/D.
 - ``owa``: outer edge of the dark zone region, in lambda_0/D.
 - ``bandwidth``: dark zone bandpass
 - ``num_wavelengths``: number of wavelengths spanning the design bandpass.
 - ``resolution``: spatial resolution of the coronagraph in each plane.

**Optimization method parameters**
 - ``force_no_x_mirror_symmetry``: force the algorithm to ignore any x mirror symmetry that might exist in the problem.
 - ``force_no_y_mirror_symmetry``: force the algorithm to ignore any y mirror symmetry that might exist in the problem.
 - ``force_no_hermitian_symmetry``: force the algorithm to ignore the (always present) Hermitian symmetry in the problem.
 - ``starting_scale``: number of pixels per unit cell for the initial solution, used in the adaptive algorithm. Must be a power of 2 times the ``ending_scale``.
 - ``ending_scale``: number of pixels per unit cell for the final solution, used in the adaptive algorithm. If this is the same as `starting_scale` the adaptive algorithm is essentially turned off.
 - ``edge_width_for_prior``: width of the optimized regions along the edges.
 - ``num_throughput_iterations``: number of iterations to let the throughput factor converge. Too few iterations and the contrast might be a fraction too high, too high and you are wasting computation time. Setting this to 1 turns iterations off.
 - ``initial_throughput_estimate``: expected relative throughput of the coronagraph compared to without apodizer or focal-plane mask (but including Lyot stop). A good estimate (for example from a lower resolution optimization) can often eliminate the need for throughput iterations altogether.
 - ``maximize_planet_throughput``:  maximize throughput through the apodizer times Lyot stop. Otherwise maximize the throughput of just the apodizer.

**Solver parameters**
 - ``num_threads``: number of threads that Gurobi is allowed to use. A value of 0 indicates that Gurobi is free to ßchoose the number.
 - ``crossover``: crossover strategy that Gurobi needs to use. See Gurobi documentation for more information.
 - ``method``: The algorithm that Gurobi needs to use. See Gurobi documentation for more information.
