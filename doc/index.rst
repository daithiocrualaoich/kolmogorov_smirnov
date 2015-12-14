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

Sequences from :math:`N(0, 1)`, :math:`N(0, 2)`, and :math:`N(1, 1)` are
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
`Github repository <https://github.com/daithiocrualaoich/kolmogorov-smirnov>`_.
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
`Github repository <https://github.com/daithiocrualaoich/kolmogorov-smirnov>`_.

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
