The Kolmogorov-Smirnov Test
===========================

.. toctree::
   :maxdepth: 2

A visit to a data and statistical technique useful to software engineers. We
learn about some Rust too along the way.

The code and examples here are available on
`Github <https://github.com/daithiocrualaoich/kolmogorov-smirnov>`_.


Datasets
--------
Statistical tests are more fun if you have datasets to run them over.

N(0,1)
~~~~~~
Because it is traditional and because it is easy and flexible, start with some
Normal distributed data.

Rust can generate Normal data using the `rand::distributions`_ module. If
``mean`` and ``variance`` are ``f64`` values representing the mean and variance
of the desired Normal deviate, then the following code generates the deviate.
Note that the ``Normal::new`` call requires the mean and standard deviation as
parameters, so it is necessary to take the square root of the variance to
provide the standard deviation value.

.. _`rand::distributions`: https://doc.rust-lang.org/rand/rand/distributions/

.. sourcecode:: rust

    extern crate rand;
    use rand::distributions::{Normal, IndependentSample};

    let mean: f64 = ...
    let variance: f64 = ...

    let mut rng = rand::thread_rng();
    let normal = Normal::new(mean, variance.sqrt);

    let x = normal.ind_sample(&mut rng);

The ``kolmogorov-smirnov`` library includes a binary for generating sequences of
independently distributed Normal deviates. It can be called with the following
usage.

::

    cargo run --bin normal <num_deviates> <mean> <variance>

The ``-q`` option is useful too in order to suppress ``cargo`` build messages in
the output.

Sequences from :math:`N(0, 1)`, :math:`N(0, 2)`, and :math:`N(1, 1)` have been
included in the
`Github repository <https://github.com/daithiocrualaoich/kolmogorov-smirnov>`_.

:math:`N(0, 2)` is largely included just to be deliberately annoying,
calculating :math:`\sqrt{2}` and drawing attention to the limitations of the
floating point represention of irrational numbers. #trololo

.. sourcecode:: bash

    cargo run -q --bin normal 8192 0 1 > normal_0_1.tsv
    cargo run -q --bin normal 8192 0 2 > normal_0_2.tsv
    cargo run -q --bin normal 8192 1 1 > normal_1_1.tsv

These are sadly not the most beautiful of Normal curves, but you must take what
you get. The :math:`N(0, 1)` data is lumpy and not single peaked.

.. image:: images/n01-1.png

:math:`N(0, 2)` is similar, though less of a disaster near the mean.

.. image:: images/n02-1.png

:math:`N(1, 1)` actally looks suprisingly like the Normal data diagrams in
textbooks.

.. image:: images/n11-1.png

And the following is a plot of all three datasets to illustrate the relative
widths, heights and supports.

.. image:: images/n01n02n11-1.png
