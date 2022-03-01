#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  cmake \
  curl \
  git \
  libbz2-dev \
  libcurl4-gnutls-dev \
  liblzma-dev \
  libssl-dev \
  mafft \
  python3-pip \
  python3-setuptools \
  tabix \
  libvcflib-tools \
  wget \
  zlib1g-dev


if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

#______________________ python ________________________________#
pip3 install tox "six>=1.14.0"

