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
RUN     pip install Sphinx
RUN     pip install cloud_sptheme


################################################################################
# Rust
################################################################################

ENV     rust_version    1.4.0-x86_64-unknown-linux-gnu
RUN     wget -q https://static.rust-lang.org/dist/rust-${rust_version}.tar.gz
RUN     tar xzf rust-${rust_version}.tar.gz
RUN     rust-${rust_version}/install.sh --without=rust-docs
RUN     rm -fr rust-${rust_version}.tar.gz rust-${rust_version}
