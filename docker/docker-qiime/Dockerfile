FROM ubuntu:14.04
# QIIME base install
MAINTAINER Yoshiki Vazquez-Baeza <yoshiki@ucsd.edu>

RUN apt-get update && apt-get -y install curl
RUN apt-get -y --no-install-recommends install cmake mercurial git make subversion
RUN apt-get -y install gcc g++ gdb
RUN apt-get -y install python python-dev python-pip vim
RUN apt-get -y install pkg-config libpng-dev libfreetype6-dev
RUN apt-get -y install libopenblas-base libopenblas-dev libatlas-base-dev gfortran libblas-dev liblapack-dev mklibs libhdf5-dev

RUN pip install numpy
RUN pip install -U qiime

ADD qiime_config /qiime_config
ENV QIIME_CONFIG_FP /qiime_config

RUN print_qiime_config.py -t

# upstream and downstream scripts to support the BaseSpace apps
ADD upstream.sh /upstream.sh
ADD downstream.sh /downstream.sh