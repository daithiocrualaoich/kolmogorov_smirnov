[![Build Status](https://travis-ci.org/daithiocrualaoich/kolmogorov-smirnov.svg?branch=master)](https://travis-ci.org/daithiocrualaoich/kolmogorov-smirnov)

Implementation of the Kolmogorov-Smirnov statistical test as a Rust library.
Read an introduction about this project, Rust, and the Kolmogorov-Smirnov test
[here](http://daithiocrualaoich.github.io/kolmogorov-smirnov).


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

To run this command directly in the Docker container without the intermediate
shell, use:

    docker run -t -v "$(pwd):/statistics" --workdir=/statistics statistics cargo test


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

    docker run -t -v "$(pwd):/statistics" --workdir=/statistics/doc statistics make clean html

The compiled document is written to the shared location and is available on the
host machine under `doc/_build`. It is published at
http://daithiocrualaoich.github.io/kolmogorov-smirnov using [Github Pages].

[Github Pages]: https://pages.github.com

To republish updated documentation, first build the html. Then create a copy of
the repository and checkout the `gh-pages` branch. A separate copy is useful
because the `master` and `gh-pages` branches are very dissimilar and switching
between them with uncommitted changes is tedious.

    cd ..
    cp -r kolmogorov-smirnov kolmogorov-smirnov-ghpages
    cd kolmogorov-smirnov-ghpages
    git reset --hard HEAD
    git clean -fdx
    git checkout gh-pages

Now remake the contents of this branch from the recently generated document.

    rm -fr *
    cp -r ../kolmogorov-smirnov/doc/_build/html/* .

A `.nojekyll` file is also needed in order to prevent Github from ignoring the
Sphinx CSS files.

    touch .nojekyll

Commit the changes and push the `gh-pages` branch to origin to perform the
publication.

    git push origin gh-pages


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
