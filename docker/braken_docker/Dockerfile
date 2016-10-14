#Creates Dockerfile for Braken
#Requires Kraken, Perl, and Python

FROM genomicpariscentre/bioperl
MAINTAINER gregpoore <gregpoore@ucsd.edu>
COPY . /usr/scripts/

# Install perl module(s) for Braken
RUN apt-get update

RUN cpanm \
 Parallel::ForkManager

# Install Python 3 via Miniconda and install useful packages

RUN apt-get update
RUN apt-get install -y xz-utils ca-certificates jq wget make g++ sudo libsm6 libxrender1 libfontconfig1 libglib2.0-0

RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.1.11-Linux-x86_64.sh --no-check-certificate && \
    bash Miniconda3-4.1.11-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm Miniconda3-4.1.11-Linux-x86_64.sh
ENV PATH /opt/miniconda3/bin:$PATH
RUN chmod -R a+rx /opt/miniconda3

# Install PyData modules and IPython dependencies
RUN conda update --yes conda && \
    conda install --yes numpy scipy pandas matplotlib cython pyzmq scikit-learn seaborn six statsmodels theano pip tornado jinja2 sphinx pygments nose readline sqlalchemy ipython jupyter


# Clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN conda clean -y -t

# Test
RUN python -c "import numpy, scipy, pandas, matplotlib, matplotlib.pyplot, sklearn, seaborn, statsmodels, theano"

# Install Kraken

ENV INSTALL_DIR /usr/bin/

# Install Kraken
RUN apt-get update
RUN apt-get install -y xz-utils ca-certificates jq wget make g++ sudo

RUN wget --output-document - http://ccb.jhu.edu/software/kraken/dl/kraken-0.10.5-beta.tgz | tar xzf - --directory $INSTALL_DIR  --strip-components=1

ENV PATH ${INSTALL_DIR}:$PATH

RUN cd $INSTALL_DIR && sudo /bin/bash install_kraken.sh $INSTALL_DIR || exit 0

RUN echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.profile
RUN echo 'export PATH="/usr/scripts:$PATH"' >> ~/.profile



