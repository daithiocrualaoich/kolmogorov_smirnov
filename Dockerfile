FROM        ubuntu:wily
MAINTAINER  Daithi O Crualaoich


################################################################################
# Basic Development Tools
################################################################################

RUN     apt-get update -qq
RUN     apt-get upgrade -qq

RUN     apt-get install -qq wget
RUN     apt-get install -qq build-essential gcc


################################################################################
# LaTeX
################################################################################

RUN     apt-get install -qq texlive
RUN     apt-get install -qq texlive-latex-extra dvipng


################################################################################
# Pandoc
################################################################################

RUN     apt-get install -qq pandoc


################################################################################
# Sphinx
################################################################################

RUN     apt-get install -qq python2.7 python2.7-dev python-pip
RUN     pip install Sphinx sphinxcontrib-googleanalytics
RUN     pip install cloud_sptheme


################################################################################
# R
################################################################################

RUN     apt-get install -qq libjpeg62 libcairo2-dev

RUN     gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN     gpg -a --export E084DAB9 | apt-key add -
RUN     echo 'deb http://cran.r-project.org/bin/linux/ubuntu wily/' > /etc/apt/sources.list.d/cran.list

RUN     apt-get install -qq r-base

RUN     echo 'update.packages(ask = FALSE, repos="http://cran.r-project.org")' | R --vanilla
RUN     echo 'install.packages(c("markdown"), repos="http://cran.r-project.org", dependencies=TRUE)' | R --vanilla
RUN     echo 'install.packages(c("knitr"), repos="http://cran.r-project.org", dependencies=TRUE)' | R --vanilla
RUN     echo 'install.packages(c("Cairo"), repos="http://cran.r-project.org", dependencies=TRUE)' | R --vanilla
RUN     echo 'install.packages(c("ggplot2"), repos="http://cran.r-project.org", dependencies=TRUE)' | R --vanilla
RUN     echo 'install.packages(c("scales"), repos="http://cran.r-project.org", dependencies=TRUE)' | R --vanilla


################################################################################
# Rust
################################################################################

ENV     rust_version    1.5.0-x86_64-unknown-linux-gnu
RUN     wget -q https://static.rust-lang.org/dist/rust-${rust_version}.tar.gz
RUN     tar xzf rust-${rust_version}.tar.gz
RUN     rust-${rust_version}/install.sh --without=rust-docs
RUN     rm -fr rust-${rust_version}.tar.gz rust-${rust_version}

RUN     cargo install --git https://github.com/rust-lang-nursery/rustfmt --root /usr
