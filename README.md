mictest
=======
    git clone git@github.com:cerati/mictest
    git checkout stable

    make -j 32

See Makefile.config for various settings. Some notes:
* icc will be detected automatically.
* to enable manually written intrinsics code (for AVX and KNC only) uncomment the following line in Section 5:
    #USE_INTRINSICS := -DMPLEX_USE_INTRINSICS
* vtune pause / resume are called to mark the core of the algorithm.

Standard test
=============
    cd mkFit
    # Generate input file
    ./mkFit --write --file-name sim-10-20k.raw --num-events 10 --num-tracks 20000
    # Run the test, e.g.
    ./mkFit --read --file-name sim-10-20k.raw --build-tbb --cloner-single-thread --num-threads 8 --best-out-of 4

With a single thread the test runs about 4 seconds. To get longer run-times either generate an input file with more events (`--num-events`) or use `--best-out-of` option to have the core of the algorithm executed several times for each event.

