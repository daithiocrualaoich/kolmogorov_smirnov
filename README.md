[![Build Status](https://travis-ci.org/daithiocrualaoich/kolmogorov_smirnov.svg?branch=master)](https://travis-ci.org/daithiocrualaoich/kolmogorov_smirnov)

Implementation of the Kolmogorov-Smirnov statistical test as a Rust library.
Read an introduction about this project, Rust, and the Kolmogorov-Smirnov test
[here](http://daithiocrualaoich.github.io/kolmogorov_smirnov).


Getting Started
---------------
The Kolmogorov-Smirnov library is available as a crate, so it is easy to
incorporate into your programs. Add the dependency to your `Cargo.toml` file.

    [dependencies]
    kolmogorov_smirnov = "1.1.0"

Information about the latest published crate is available on
[crates.io](https://crates.io/crates/kolmogorov_smirnov).

Using the test is also straightforward, call the `kolmogorov_smirnov::test`
function with the two samples to compare and the desired confidence level.

    extern crate kolmogorov_smirnov as ks;

    let xs = vec!(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
    let ys = vec!(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    let confidence = 0.95;

    let result = ks::test(&xs, &ys, confidence);

    if !result.is_rejected {
        // Woot! Samples are from the same distribution with 0.95 confidence.
    }

Alternatively, if you have floating point or integer data to test, you can use
the included test runners, ``ks_f64.rs`` and ``ks_i32.rs``. These operate on
single-column headerless data files and test the samples against each other at
the 0.95 confidence level.

    $ cargo run -q --bin ks_f64 dat/normal_0_1.tsv dat/normal_0_1.1.tsv
    Samples are from the same distribution.
    test statistic = 0.0399169921875
    critical value = 0.08550809323787689
    reject probability = 0.18365715210599798

    $ cargo run -q --bin ks_f64 dat/normal_0_1.tsv dat/normal_1_1.1.tsv
    Samples are from different distributions.
    test statistic = 0.361572265625
    critical value = 0.08550809323787689
    reject probability = 1


Developing Kolmogorov-Smirnov
-----------------------------
Install the [Rust] development tools on your system if they are not already
available. Then build and test the library using:

    cargo test

[Rust]: https://www.rust-lang.org


Docker
------
A [Docker] container definition is provided with installations of the tools
used to develop the software. To use the container, first install Docker if not
already available and start a Docker terminal. Then create the container by
running the following build at the top level of the repository source tree:

    docker build --rm=true -t statistics .

[Docker]: http://docker.io

Once built, an interactive shell can be run in the container using:

    docker run -it -v "$(pwd):/statistics" --workdir=/statistics statistics /bin/bash

The current working directory from the host machine is available as the current
directory in the container so it is possible to build and test the library as
described earlier.

    cargo test

To see more detail on panics include the `RUST_BACKTRACE` environment variable.

    RUST_BACKTRACE=1 cargo test


Building the Documentation
--------------------------
The RestructuredText format [Sphinx] documentation under `doc` can be compiled
using the Makefile.

    cd doc
    make clean html

[Sphinx]: http://sphinx-doc.org

See this [RestructuredText Primer] for guidance on writing RestructuredText.

[RestructuredText Primer]: http://sphinx-doc.org/rest.html

The Docker container provides an installation of Python, Sphinx, and LaTeX
required to do this build. To make the documentation directly in container
without an intermediate shell, use:

    docker run -v "$(pwd):/statistics" --workdir=/statistics/doc statistics make clean html

The compiled document is written to the shared location and is available on the
host machine under `doc/_build`. It is published at
http://daithiocrualaoich.github.io/kolmogorov_smirnov using [Github Pages].

[Github Pages]: https://pages.github.com

To republish updated documentation, first build the html. Then create a copy of
the repository and checkout the `gh-pages` branch. A separate copy is useful
because the `master` and `gh-pages` branches are very dissimilar and switching
between them with uncommitted changes is tedious.

    cd ..
    cp -r kolmogorov_smirnov kolmogorov_smirnov_ghpages
    cd kolmogorov_smirnov_ghpages
    git reset --hard HEAD
    git clean -fdx
    git checkout gh-pages

Now remake the contents of this branch from the recently generated document.

    rm -fr *
    cp -r ../kolmogorov_smirnov/doc/_build/html/* .

A `.nojekyll` file is also needed in order to prevent Github from ignoring the
Sphinx CSS files.

    touch .nojekyll

Commit the changes and push the `gh-pages` branch to origin to perform the
publication.

    git push origin gh-pages

To regenerate the diagrams, run the following.

    docker run -v "$(pwd):/statistics" --workdir=/statistics/dat statistics R -e "rmarkdown::render('images.Rmd')"


Publishing on crates.io
-----------------------
Instructions for uploading to the crate repository at crates.io are
[here](http://doc.crates.io/crates-io.html#publishing-crates). First login to
the site using:

    cargo login <token>

Token can be found from [crates.io/me](https://crates.io/me). To make a release,
first clean and build the package:

    git stash
    cargo clean
    cargo package

Examine the built package under `target/package/kolmogorov_smirnov-<version>`.
And when happy to publish:

    cargo publish

And check out the new update at
[crates.io/crates/kolmogorov_smirnov](https://crates.io/crates/kolmogorov_smirnov).


License
-------

    Copyright [2015] [Daithi O Crualaoich]

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
