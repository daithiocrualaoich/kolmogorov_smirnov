The Kolmogorov-Smirnov Test
===========================

.. toctree::
   :maxdepth: 2

A visit to a data and statistical technique useful to software engineers. We
learn about some Rust too along the way.

The code and examples here are available on
`Github <https://github.com/daithiocrualaoich/kolmogorov_smirnov>`_.


Kolmogorov-Smirnov Hypothesis Testing
-------------------------------------
The Kolmogorov-Smirnov test is a hypothesis test procedure for determining if
two samples of data come from the same distribution. The test is non-parametric,
entirely agnostic to what this distribution actually is. The fact that we never
have to know what distribution the samples come from is incredibly useful,
especially in software and operations where the distributions are hard to
express and difficult to calculate with.

It is really surprising that such a test exists. In an unkind Universe, we would
be completely on our own.

The test description may look a bit hard in the outline below. But skip ahead to
the implementation because the Kolmogorov-Smirnov test is incredibly easy in
practice.

The Kolmogorov-Smirnov test is covered in `Numerical Recipes`_. There is a
`pdf <http://www.aip.de/groups/soe/local/numres/bookcpdf/c14-3.pdf>`_ available
from the third edition of Numerical Recipes in C.

.. _Numerical Recipes: http://www.aip.de/groups/soe/local/numres

The `Wikipedia article <https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test>`_
is a useful overview but light about proof details. If you are interested in a
why the test statistic has a distribution that is independent and useful for
constructing the test, then these
`MIT lecture notes <http://ocw.mit.edu/courses/mathematics/18-443-statistics-for-applications-fall-2006/lecture-notes/lecture14.pdf>`_
give a sketch overview.

For an application to metrics and monitoring in software operations, see this
`introductory talk <https://vimeo.com/95069158>`_ by Toufic Boubez at
Monitorama. The slides are available on
`slideshare <http://www.slideshare.net/tboubez/simple-math-for-anomaly-detection-toufic-boubez-metafor-software-monitorama-pdx-20140505>`_.

The Test Statistic
~~~~~~~~~~~~~~~~~~
The Kolmogorov-Smirnov test is constructed as a statistical hypothesis test. We
determine a null hypothesis, :math:`H_0`, that the two samples we are testing
come from the same distribution. Then we search for evidence that this
hypothesis should be rejected and express this in terms of a probability. If
the likelihood of the samples being from different distributions exceeds a
confidence level we demand then the original hypothesis is rejected in favour
of the hypothesis, :math:`H_1`, that the two samples are from different
distributions.

To do this, we devise a single number calculated from the samples, i.e. a
statistic. The trick is to find a number that has a range of values which do
not depend on things we don't know like the actual underlying distributions in
this case.

The test statistic in the Kolmogorov-Smirnov test is very easy, just the
maximum vertical distance between the empirical cumulative distribution
functions of the two samples.

For instance, in this plot of the empirical cumulative distribution functions
of a Normal(0, 1) and a Normal(0, 2) sample, the maximum vertical distance
between the lines is at about -1.5 and 1.5.

.. image:: images/n01n02ecdf-1.png

For a Normal(0, 1) against Normal(1, 1) sample, it is a lot clearer. The maximum
vertical distance occurs somewhere around zero and is quite large, maybe about
0.35 in size.

.. image:: images/n01n11ecdf-1.png

As an aside, these demonstrate an important note about the application of the
Kolomogorov-Smirnov test. It is much better at detecting distributions where
the medians are far apart than it is at detecting distributions where the tails
are different but the main mass of the distributions are around the same values.

So, if :math:`X_i` are n independent and identically distributed observations,
the empirical cumulative distribution function :math:`F_n` is:

.. math::

    F_n(x) = \frac{1}{n}\sum_{i=1}^n I_{[-\infty,x]}(X_i)

Here, :math:`I` is the indicator function which is 1 if :math:`X_i` is less than
or equal to :math:`x` and 0 otherwise.

This just says that :math:`F_n(x)` is the number of samples observed that are
less than or equal to :math:`x` divided by the total number of samples. But it
does it in a complicated way so we can feel clever about ourselves.

The empirical cumulative distribution function is an unbiased estimator for the
underlying cumulative distribution function, incidentally.

For two samples having empirical cumulative distribution functions
:math:`F_n(x)` and :math:`G_m(x)`, the Kolmogorov-Smirnov test statistic,
:math:`D`, is the maximum absolute difference between these for the same
:math:`x`,  i.e. the largest vertical distance between the plots in the graph.

.. math::

    D = \sup_{-\infty\leq x \leq\infty} |F_n(x) - G_m(x)|

The Glivenko–Cantelli theorem says if the :math:`F_n(x)` is made from samples
from the same distribution as :math:`G_m(x)` then this statistic "almost surely
converges to 0 in the limit when n goes to infinity." This is an extremely
technical statement that we are going to ignore.

Two-Sample KS test
~~~~~~~~~~~~~~~~~~
Surprisingly, the distribution of :math:`D` can be approximated well in the case
that the samples are drawn from the same distribution. This means we can build
a statistic test that rejects this null hypothesis for a given confidence level
if :math:`D` exceeds an easily calculable value.

Tables of critical values are available. The table
`here <https://www.webdepot.umontreal.ca/Usagers/angers/MonDepotPublic/STT3500H10/Critical_KS.pdf>`_
describes a test implementation for samples sizes greater than twelve where we
reject the null hypothesis, i.e. decide that the samples are from different
distributions, if

.. math::

    D > c(\alpha)\sqrt{\frac{n + m}{n m}}

Where n and m are the sample sizes. For :math:`\alpha = 0.05` corresponding to
95% confidence level, the value is :math:`c(\alpha) = 1.36`.

Numerical Recipes describes a direct calculation that works well for:

.. math::

    N_{n, m} = \frac{n + m}{n m} \geq 4

i.e. for sample sizes greater than seven as :math:`N_{8, 8} = 4`.

The probability that the test statistic :math:`D` is greater than the value
observed is approximately:

.. math::

    P(D > \text{observed}) = Q_{KS}\Big(\Big[\sqrt{N_{n, m}} + 0.12 + 0.11/\sqrt{N_{n, m}}\Big] D\Big)

where

.. math::

    Q_{KS}(x) = 2 \sum_{j=1}^{\infty} (-1)^{j-1} e^{-2j^2x^2}

The implementation in Numerical Recipes gives this a hundred terms to converge
before failing.

Discussion
~~~~~~~~~~
The implementation of this test is straightforward and can be found in the
`Github repository <https://github.com/daithiocrualaoich/kolmogorov_smirnov>`_.
The empirical cumulative distribution function implementation is about as
complicated as it gets for this. There are two functions for this in the code,
the simpler version used to probabilistically verify the other.

The non-parametricity and generality is the great advantage of the
Kolomogorov-Smirnov test but this is balanced by a number of drawbacks in its
ability to establish evidence to reject the null hypothesis.

In particular, the Kolmogorov-Smirnov test is weak in the cases that the sample
empirical distribution functions do not deviate strongly even though the samples
are from different distributions. For instance, the Kolomogorov-Smirnov test
is most sensitive around the median of the samples because this is where
differences in the graph are most likely to be. It is less strong near the tails
because both cumulative distribution functions will be near 0 or 1 and the
differences between them less pronounced. Location and shape related scenarios
that limit the :math:`D` test statistic reduce the ability of the
Kolmogorov-Smirnov test to reject the null hypothesis.

The Chi-squared test is also used for testing if samples come from the same
distribution but this is done with a binning and discretizes the data. This is
not a issue in the Kolomogorov-Smirnov test.


A Field Manual for Rust
-----------------------
`Rust`_ is a Mozilla sponsored project to create a safe, fast systems language.

.. _Rust: https://www.rust-lang.org

Why? There is a entire free
`O'Reilly book <http://www.oreilly.com/programming/free/files/why-rust.pdf>`_
on exactly this question but the reasons include:

* Robust memory management. It is not possible to deference null or dangling
  pointers in Rust.
* Improved security, reducing the incidence of flaws like buffer overflow
  exploits.
* A light runtime, no garbage collection, and overhead means Rust is convenient
  to embed in other languages and platforms like Ruby, Python, and Node.
* Modern language features.

Rust is capable of very serious projects. The current flagship Rust project, for
instance, is `Servo`_, a browser engine under open source development with
contributions from Mozilla and Samsung.

.. _Servo: https://servo.org

The best introduction to Rust is the `Rust Book`_. But also read Steve Klabnik's
`alternative introdution to Rust`_ for the upfront no-nonsense dive into memory
ownership, the crux concept for Rust beginners.

.. _Rust Book: https://doc.rust-lang.org/book
.. _alternative introdution to Rust: http://words.steveklabnik.com/a-new-introduction-to-rust

Those in a hurry can quickstart with these slide decks:

* `Dimiter Petrov, Romain Ruetschi <http://hackepfl.github.io/rust-workshop-slides/slides.html>`_
* `Danilo Bargen <https://speakerdeck.com/dbrgn/intro-to-rust>`_

Two must-read learning resources are `24 Days of Rust`_, a charming tour around
the libraries and world of Rust, and `ArcadeRS`_, a tutorial via the medium of
video games.

.. _24 Days of Rust: http://zsiciarz.github.io/24daysofrust
.. _ArcadeRS: http://jadpole.github.io/arcaders/arcaders-1-0

And finally, if Servo has interested you in the idea of writing a browser engine
in Rust, then `Let's build a browser engine!`_ is the series for you. It walks
through writing a simple HTML rendering engine in Rust.

.. _Let's build a browser engine!: http://limpet.net/mbrubeck/2014/08/08/toy-layout-engine-1.html

Moral Support for Learning the Memory Rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There is no pretending otherwise, the Road to Rust is not royal. The Rust memory
rules about lifetime, ownership, and borrowing are especially hard to learn.

It probably doesn't much feel like it but Rust is really trying to help us with
these rules. And to be fair to Rust, it hasn't segfaulted me so far.

But that is not much comfort when the compiler won't build your code and you
have no idea why. The best advice is probably to read as much about the Rust
memory rules as you can and to keep reading about them over and over until they
start to make some sense.

Adherence to the rules provides the compiler with invariant guarantees that can
be used to construct a proof of memory safety. In practical terms, though, the
rationale for these rules is unimportant. What is necessary is to find a way to
work with them so that your programmes compile.

Keep at it. It takes a long time but eventually it does all become clearer.

Niche Observations
~~~~~~~~~~~~~~~~~~
This section is a random sample of Rust arcana that interested me. Nothing here
that doesn't interest you is worth troubling too much with and you should skip
on past.

`Travis CI`_ has excellent support for building Rust projects, including with
the beta, and nightly versions. It is simple to set up too by configuring a
``travis.yml`` according to the
`Travis documentation <https://docs.travis-ci.com/user/languages/rust>`_. The
Travis CI build for this project, for instance, is
`here <https://travis-ci.org/daithiocrualaoich/kolmogorov_smirnov>`_.

.. _Travis CI: https://travis-ci.org

Rust has a formatter in `rustfmt`_ and a lint in `rust-clippy`_. The formatter
is a simple install using ``cargo install`` and provides a binary command. The
lint requires more integration into your project, and currently also needs the
nightly version of Rust for plugin support. Both projects are great for helping
Rust newcomers.

.. _rustfmt: https://github.com/rust-lang-nursery/rustfmt
.. _rust-clippy: https://github.com/Manishearth/rust-clippy

Foreign Function Interface is an area where Rust excels. The absence of a large
runtime means Rust is great for embedding in other languages and it has a wide
range as a C replacement in writing modules for Python, Ruby, Node, etc. See the
`Rust Book introduction <https://doc.rust-lang.org/book/rust-inside-other-languages.html>`_
for a demonstration of how easy it is call Rust from other languages.
`Day 23 of Rust <https://zsiciarz.github.io/24daysofrust/book/day23.html>`_ and
the `Rust FFI Omnibus`_ are additional resources for Rust FFI.

.. _Rust FFI Omnibus: http://jakegoulding.com/rust-ffi-omnibus

Rust is being used experimentally for embedded development. `Zinc`_ is work on
building a realtime operating system for ARM using Rust primarily, and the
following are posts about building software for embedded devices directly using
Rust.

* `Rust bare metal on ARM microcontroller`_.
* `Embedded Rust Right Now!`_

.. _Zinc: http://zinc.rs
.. _Rust bare metal on ARM microcontroller: http://antoinealb.net/programming/2015/05/01/rust-on-arm-microcontroller.html
.. _Embedded Rust Right Now!: http://spin.atomicobject.com/2015/02/20/rust-language-c-embedded

And relatedly, `Rust on Raspberry Pi`_ is a guide to cross-compiling Rust code
for the Raspberry Pi.

.. _Rust on Raspberry Pi: https://github.com/Ogeon/rust-on-raspberry-pi/

Rust considers the code snippets in your project documentation as tests and
makes a point of compiling them. This helps a little to keep your documentation
in sync with the code but it is a shock the first time you get a compiler error
for a documentation code snippet and it takes you ages to find it.


Kolmogorov-Smirnov Test Library
-------------------------------
The Kolmogorov-Smirnov test implementation is available as a crate, so it is
easy to incorporate into your programs. Add the dependency to your
``Cargo.toml`` file.

::

    [dependencies]
    kolmogorov_smirnov = "0.1.0"

Then using the test is straightforward, call the ``kolmogorov_smirnov::test``
function with the two samples to compare and the desired confidence level.

.. sourcecode:: rust

    extern crate kolmogorov_smirnov as ks;

    let xs = vec!(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
    let ys = vec!(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    let confidence = 0.95;

    let result = ks::test(&xs, &ys, confidence);

    if !result.is_rejected {
        // Woot! Samples are from the same distribution with 0.95 confidence.
    }

The Kolmogorov-Smirnov test as implemented works for any data with a ``Clone``
and an ``Ord`` trait implementation in Rust. It would be possible, maybe a
little eccentric, to test two samples of characters, strings or lists.

If you have floating point or integer data to test, you can use the included
test runners, ``ks_f64.rs`` and ``ks_i32.rs``. These operate on single-column
headerless data files and test the samples against each other at the 0.95
confidence level.

.. sourcecode:: bash

    $ cargo run -q --bin ks_f64 dat/normal_0_1.tsv dat/normal_0_1.1.tsv
    Samples are from the same distributions.

    $ cargo run -q --bin ks_f64 dat/normal_0_1.tsv dat/normal_1_1.1.tsv
    Samples are from different distributions.

Testing floating point numbers is a real headache because Rust floating point
types do not implement the ``Ord`` trait, only the ``PartialOrd`` trait. This
is because things like ``NaN`` are not comparable and the order cannot be
total over all values in the datatype.

The test runner for floating point types is implemented with a wrapper type
that implements a total order, crashing on unorderable elements. This suffices
in practice since the unorderable elements will break the test anyway.


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

The ``kolmogorov_smirnov`` library includes a binary for generating sequences of
independently distributed Normal deviates. It can be called with the following
usage.

::

    cargo run --bin normal <num_deviates> <mean> <variance>

The ``-q`` option is useful too in order to suppress ``cargo`` build messages in
the output.

Sequences from :math:`N(0, 1)`, :math:`N(0, 2)`, and :math:`N(1, 1)` are
included in the
`Github repository <https://github.com/daithiocrualaoich/kolmogorov_smirnov>`_.

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

http://twitter.com Response Times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Less artificially, and more likely in software operations and monitoring, are
datasets of HTTP server response times. Metrics behaving like response time are
very typical and not easy to analyse with standard statistical technology.

`Apache Bench`_ is a commandline URL loadtesting tool that can be used to
collect sample datasets for HTTP request times. A dataset of service times for
``http://twitter.com`` was collected using:

.. _Apache Bench: https://httpd.apache.org/docs/2.2/programs/ab.html

.. sourcecode:: bash

    ab -n 8192 -c 1 -v 3 -g http.tsv http://twitter.com/

The test actually ships 3.36MB of redirect headers since ``http://twitter.com``
is a 301 Moved Permanently redirect to ``https://twitter.com``. The dataset is
still useful as it exhibits the behaviour of HTTP endpoint responses anyway.

(Note that Apache Bench has an option for issuing HEAD requests instead of GET
requests so it is possible to test the HTTPS endpoint without the 2GB of HTML
traffic but the encryption setup still has a cost for everybody.)

The options in the call specify:

* ``-c 1``: Use test concurrency of one outstanding request at any given time to
  trottle the testing.
* ``-v 3``: Log at a high verbosity level to show the individual request data.
* ``-g http.tsv``: Output a TSV summary to ``http.tsv``. The ``-g`` stands for
  Gnuplot which the output file is specifically organised to support.

A timeout or significant failure can trash the test completely, so it is more
robust to collect the data in blocks of 256 requests and combine the results.
This was done to collect some supplementary data for comparison purposes and is
available as ``dat/http.1.tsv`` through ``dat/http.4.tsv`` in the
`Github repository <https://github.com/daithiocrualaoich/kolmogorov_smirnov>`_.
The primary dataset is ``dat/http.tsv``.

.. sourcecode:: bash

    for i in {1..32}
    do
      ab -n 256 -c 1 -v 3 -g http.$i.tsv http://twitter.com/
    done

Aside: The trailing ``/`` in the target URL is required in the ``ab`` command or
it fails with an invalid URL error.

The result data file from Apache Bench has the following
`schema <http://stackoverflow.com/a/6143973>`_:

* ``starttime``: A human friendly time representation of the request start.
* ``seconds``: The Unix epoch seconds timestamp of ``starttime``. This is the
  number of regular, but not leap, seconds since 1 January 1970. Trivia Fact:
  this is why Unix timestamps line up, against sane expectations, on minute
  boundaries when divided by 60.
* ``ctime``: Connection time to the server in milliseconds.
* ``dtime``: Processing time on the server in milliseconds. The ``d`` may stand
  for "duration" or it may not.
* ``ttime``: Total time in milliseconds, ``ctime + dtime``
* ``wait``: Waiting time in milliseconds. This is not included in ``ttime``, it
  appears from the data.

If you want to the generated data file as it was intended for processing with
Gnuplot, then see this
`article <http://www.bradlanders.com/2013/04/15/apache-bench-and-gnuplot-youre-probably-doing-it-wrong>`_.

The output data is sorted by Apache Bench according to ``ttime``, so shorter
requests come first in the output. The purposes here better appreciate data
sorted by the ``seconds`` timestamp, particularly to plot the timeseries data in
the order it happened.

And because nobody is much interested in when exactly this dataset was made, the
``starttime`` and ``seconds`` fields can be dropped from the final data. The
following `crafty piece of Awk <http://unix.stackexchange.com/a/71949>`_ does
the desired header-retaining sort and column projection.

.. sourcecode:: bash

    cat http.tsv |\
      awk 'NR == 1; NR > 1 {print $0 | "sort"}' |\
      cut -f3- \
      > http.tsv.out

The response time datasets in the Github repository have all been processed like
this.

To process the headered multi-column data files into a format suitable for
using with ``ks_i64`` use the following example as a guide:

.. sourcecode:: bash

    tail -n +2 http.tsv |\
      cut -f3 \
      > http_ttime.tsv

Then to run the test:

.. sourcecode::bash

    cargo run -q --bin test_t64 dat/http_ttime.tsv dat/http_ttime.1.tsv

The timeseries plot shows a common Internet tale of outliers, failure cases and
disgrace.

.. image:: images/http-timeseries-1.png

The density is highly peaked but includes fast failure weight and a long light
tail. This is not straightforward to parametrise with common statistical
distributions that are fruitful to work with and demonstrates a significant
utility of the Kolmogorov-Smirnov test.

.. image:: images/http-density-1.png

Note there is a much larger horizontal axis range in this graph. This has the
effect of compressing the visual area under the graph relative to the other
dataset density plots.

Don't let this trick you, though. The y axis values are smaller than in the
other graphs but there is far more horizontal axis support to compensate. The
definition of a probability density means that the area under the graph in all
the density plots must sum to 1.

Twitter, Inc. Stock Price
~~~~~~~~~~~~~~~~~~~~~~~~~
The final dataset is a historical stock market sample. The following returns a
14 day dump of minute granularity stock price data for Twitter, Inc. from Google
Finance.

.. sourcecode:: bash

    wget -O twtr.dat 'http://www.google.com/finance/getprices?i=60&p=14d&f=d,c,h,l,o,v&q=TWTR'

The options in the call specify:

* `i=60`: Sample interval in seconds, i.e. get per-minute data. The minimum
  sample interval available is 60 seconds.
* `p=14d`: Show data for the previous 14 days.
* `f=d,c,h,l,o,v`: Include columns in the result for sample interval start date,
  sample interval closing price, sample interval high price value, sample
  interval low price value, sample interval opening price, and sample interval
  trade count, i.e. volume.
* `q=TWTR`: Query data for the TWTR stock symbol.

After the header block in the response, the sample data looks like:

::

    a1448289000,26.11,26.11,26.11,26.11,500
    1,26.17,26.17,26.12,26.13,700
    2,26.25,26.25,26.14,26.14,191266
    3,26.14,26.27,26.13,26.21,89148
    4,26.01,26.15,26.01,26.135,36535

Lines that start with an ``a`` character include an absolute timestamp value.
Otherwise, the timestamp field value is an offset and must be added to timestamp
value in the last previous absolute timestamp line.

The following Awk script truncates the header block and folds the timestamp
offsets into absolute timestamp values.

.. sourcecode:: bash

    tail +8 twtr.dat | awk '
      BEGIN {
        FS=","; OFS="\t"
        print "timestamp","close","high","low","open","volume"
      } {
        if (substr($1,0,1) == "a") {
          $1 = substr($1,2)
          base = 0 + $1
        } else {
          $1 = base + (60 * $1)
        }

        print $0
      }' > twtr.tsv

The output TSV file is available as ``dat/twtr.tsv`` in the
`Github repository <https://github.com/daithiocrualaoich/kolmogorov_smirnov>`_.

A supplementary dataset consisting of a single day was collected for comparison
purposes and is available as ``dat/twtr.1.tsv``. It was processed identically to
``dat/twtr.tsv``. The collection command was:

.. sourcecode:: bash

    wget -O twtr.1.dat 'http://www.google.com/finance/getprices?i=60&p=1d&f=d,c,h,l,o,v&q=TWTR'

Trading hours on the New York Stock Exchange are weekdays 9.30am to 4pm. This
results in long horizontal line segments for missing values in the timeseries
plot for interval opening prices. These correspond to the overnight and weekend
market close periods.

.. image:: images/twtr-open-timeseries-1.png

The following is the opening price timeseries density plot.

.. image:: images/twtr-open-density-1.png

The missing value graph artifact is more pronounced in the trading volume
timeseries. The lines joining the trading day regions can be disregarded.

.. image:: images/twtr-volume-timeseries-1.png

Finally, the trading volume density plot looks surprisingly structured,
congregating near the 17,000 trades/minute rate.

.. image:: images/twtr-volume-density-1.png


A Diversion In QuickCheck
-------------------------
`QuickCheck`_ is insane amounts of fun writing tests and a great way to become
comfortable in a new language.

.. _QuickCheck: https://en.wikipedia.org/wiki/QuickCheck

The QuickCheck idea is to write tests as properties of the inputs rather than
specific test cases. So, for instance, rather than checking whether a given
pair of samples have a determined maximum empirical cumulative distribution
function distance, instead a generic property is verified such as the distance
is between zero and one for any pair of input samples.

This form of test construction means that QuickCheck can probabilistically
check the property over a huge number of test case instances and establish much
greater confidence of correctness than a single individual test instance could.

It can be harder too, yes. Writing properties that tightly specify the desired
behaviour is an art form but starting with properties that very loosely
constrain the software behaviour is often helpful, and facilitates an evolution
into more sharply binding criteria.

For a tutorial introduction to QuickCheck, John Hughes has a great
`introduction talk <https://www.youtube.com/watch?v=zi0rHwfiX1Q>`_.

There is an implementation of `QuickCheck for Rust`_ and the tests for the
Kolmogorov-Smirnov Rust library have been implemented using it. See the
`Github repository <https://github.com/daithiocrualaoich/kolmogorov_smirnov>`_
for examples of how to QuickCheck in Rust.

.. _QuickCheck for Rust: https://github.com/BurntSushi/quickcheck

Here is a toy example of a QuickCheck property to test an integer doubling
function.

.. sourcecode:: rust

    extern crate quickcheck;
    use self::quickcheck::quickcheck;

    fn double(n: u32) -> u32 {
        2 * n
    }

    #[test]
    fn test_double_n_is_greater_than_n() {
        fn prop(n: u32) -> bool {
            double(n) > n
        }

        quickcheck(prop as fn(u32) -> bool);
    }

This test is broken and QuickCheck makes short(I almost wrote 'quick'!) work of
letting us know that we have been silly.

::

    test tests::test_double_n_is_greater_than_n ... FAILED

    failures:

    ---- tests::test_double_n_is_greater_than_n stdout ----
        thread 'tests::test_double_n_is_greater_than_n' panicked at '[quickcheck] TEST FAILED. Arguments: (0)', /root/.cargo/registry/src/github.com-0a35038f75765ae4/quickcheck-0.2.24/src/tester.rs:113

The last log line here includes the ``u32`` instance that failed the test, zero.
Correct practice is to now create a non-probabilistic test case from the failure
that tests that specific value. This will protect the codebase from regressions
in the future.

The problem is that the property is not actually valid for the ``double``
function because double zero is actually zero. So let's fix it.

.. sourcecode:: rust

    #[test]
    fn test_double_n_is_geq_n() {
        fn prop(n: u32) -> bool {
            double(n) >= n
        }

        quickcheck(prop as fn(u32) -> bool);
    }

Note how QuickCheck produced a minimal test violation, there are no smaller
values of ``u32`` that violated the test. This is not an accident. QuickCheck
libraries often include support for shrinking test failures to minimal examples.
So when a test fails, QuickCheck will rerun it searching successively on smaller
instances of the test arguments to determine the smallest violating case.

The function is still broken, by the way, because it will overflow for large
integers. QuickCheck doesn't catch this problem because the
``QuickCheck::quickcheck`` convenience runner configures the tester to produce
random data between zero and one hundred. For this reason, you should not use
this convenience runner in testing. Instead, configure the ``QuickCheck``
manually with as large a random range as you can.

.. sourcecode:: rust

    extern crate rand;
    use self::quickcheck::{QuickCheck, StdGen, Testable};
    use std::usize;

    fn check<A: Testable>(f: A) {
        let gen = StdGen::new(rand::thread_rng(), usize::MAX);
        QuickCheck::new().gen(gen).quickcheck(f);
    }

    #[test]
    fn test_double_n_is_greater_than_n_if_n_is_greater_than_1() {
        fn prop(n: u32) -> TestResult {
            if n <= 1 {
                return TestResult::discard()
            }

            TestResult::from_bool(double(n) > n)
        }

        check(prop as fn(u32) -> TestResult);
    }

This will break the test with an overflow panic. This is correct and the
``double`` function should do something about handling overflow properly.

A warning, though, if you are testing ``vec``s or strings. The number of
elements in the randomly generated ``vec`` or equivalently, the length of the
randomly generated string will be a random number between zero and the size in
the ``StdGen`` configured. There is the potential in this to create
unnecessarily huge ``vec``s and strings. See the example of ``NonEmptyVec``
below for a technique to limit the size of a randomly generated ``vec`` or
string while still using a ``StdGen`` with a large range.

Unfortunately, you are out of luck on a 32-bit machine where the ``usize::MAX``
will only get you to sampling correctly in ``u32``. You will need to upgrade to
a new machine before you can test ``u64``, sorry.

By way of example, it is actually more convenient to include known failure cases
like ``u32::max_value()`` in the QuickCheck test function rather than in a
separate traditional test case function because the property code is available.
So, when the QuickCheck fails for the overflow bug, add the test case like
follows instead of in a separate function:

.. sourcecode:: rust

    #[test]
    fn test_double_n_is_geq_n() {
        fn prop(n: u32) -> bool {
            double(n) >= n
        }

        assert!(prop(u32::max_value()));
        quickcheck(prop as fn(u32) -> bool);
    }

Sometimes the property to test is not valid for some test arguments, i.e. the
property is useful to verify but there certain combinations of probabilistically
generated inputs that should be excluded.

The Rust QuickCheck library supports this with the ``TestResult`` type. Suppose
that instead of writing the ``double`` test property correctly, we wanted to
just exclude the failing cases instead. This might be a practical thing to do in
a real scenario. To do this, we rewrite the test as follows:

.. sourcecode:: rust

    use self::quickcheck::TestResult;

    #[test]
    fn test_double_n_is_greater_than_n_if_n_is_greater_than_1() {
        fn prop(n: u32) -> TestResult {
            if n <= 1 {
                return TestResult::discard()
            }

            TestResult::from_bool(double(n) > n)
        }

        quickcheck(prop as fn(u32) -> TestResult);
    }

Here the cases where the property legimately doesn't hold are excluded by
returning a ``TestResult::discard()``. This causes QuickCheck to retry the test
with the next randomly generated ``u32``.

Note also that the function return type is now ``TestResult`` and that
``TestResult::from_bool`` is required around the test condition.

An alternative approach is to create a wrapper type in the test code which only
permits valid input and to rewrite the tests to take this type as the
probabilistically generated input instead.

For example, suppose you want to ensure that QuickCheck only generates positive
integers for use in your property verification, you could add a wrapper type
``PositiveInteger``. For QuickCheck to work, you then have to implement the
``Arbitrary`` trait for this new type.

The minimum requirement for an ``Arbitrary`` implementation is a function called
``arbitrary`` taking a ``Gen`` random generator and producing a
``PositiveInteger``. New implementations should always leverage existing
``Arbitrary`` implementations, and so ``PositiveInteger`` generates a random
``u64`` using ``u64::arbitrary()`` and constrains it to be greater than zero.

.. sourcecode:: rust

    extern crate quickcheck;
    use self::quickcheck::{Arbitrary, Gen};
    use std::cmp;

    #[derive(Clone, Debug)]
    struct PositiveInteger {
        val: u64,
    }

    impl Arbitrary for PositiveInteger {
        fn arbitrary<G: Gen>(g: &mut G) -> PositiveInteger {
            let val = cmp::max(u64::arbitrary(g), 1);
            PositiveInteger { value: val }
        }

        fn shrink(&self) -> Box<Iterator<Item = PositiveInteger>> {
            let shrunk: Box<Iterator<Item = u64>> = self.value.shrink();
            Box::new(shrunk.filter(|&v| v > 0).map(|v| {
                PositiveInteger { value: v }
            }))
        }
    }

Note also the implementation of ``shrink()`` here, again in terms of the
existing ``u64::shrink()``. This method is optional, if it is unimplemented then
QuickCheck won't minimise property violations for the new wrapper type.

Use ``PositiveInteger`` like follows:

.. sourcecode:: rust

    use self::quickcheck::quickcheck;

    fn square(n: u64) -> u64 {
        n * n
    }

    #[test]
    fn test_square_n_for_positive_is_geq_1() {
        fn prop(n: PositiveInteger) -> bool {
            square(n.value) >= 1
        }

        quickcheck(prop as fn(PositiveInteger) -> bool);
    }

There is no need to use ``TestResult::discard()`` to ignore the failure case
where n is zero.

Finally, wrappers can be added for more complicated types too. A commonly
useful container type generator is ``NonEmptyVec`` which excludes the empty
``vec``.

.. sourcecode:: rust

    extern crate quickcheck;
    extern crate rand;
    use self::quickcheck::{quickcheck, Arbitrary, Gen};
    use self::rand::Rng;
    use std::cmp;

    #[derive(Debug, Clone)]
    struct NonEmptyVec<A> {
        value: Vec<A>,
    }

    impl<A: Arbitrary> Arbitrary for NonEmptyVec<A> {
        fn arbitrary<G: Gen>(g: &mut G) -> NonEmptyVec<A> {
            // Limit size of generated vec to 1024
            let max = cmp::min(g.size(), 1024);

            let size = g.gen_range(1, max);
            let vec = (0..size).map(|_| A::arbitrary(g)).collect();

            NonEmptyVec { value: vec }
        }

        fn shrink(&self) -> Box<Iterator<Item = NonEmptyVec<A>>> {
            let vec: Vec<A> = self.value.clone();
            let shrunk: Box<Iterator<Item = Vec<A>>> = vec.shrink();

            Box::new(shrunk.filter(|v| v.len() > 0).map(|v| {
                NonEmptyVec { value: v }
            }))
        }
    }

    #[test]
    fn test_head_of_sorted_vec_is_smallest() {
        fn prop(vec: NonEmptyVec<u64>) -> bool {
            let mut sorted = vec.value.clone();
            sorted.sort();

            // NonEmptyVec must have an element.
            let head = sorted[0];

            vec.value.iter().all(|&n| head <= n)
        }

        quickcheck(prop as fn(NonEmptyVec<u64>) -> bool);
    }
