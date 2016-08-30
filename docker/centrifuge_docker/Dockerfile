# Install Centrifuge

FROM ubuntu:14.04

ENV INSTALL_DIR /usr/local/bin/

RUN apt-get update
RUN apt-get install -y unzip xz-utils ca-certificates jq wget make g++ sudo bsdtar

RUN wget --output-document - ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.2-beta-Linux_x86_64.zip | bsdtar xvf - -C $INSTALL_DIR --strip-components=1

ENV PATH ${INSTALL_DIR}:$PATH

RUN cd $INSTALL_DIR && make $INSTALL_DIR || exit 0

RUN echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.profile

RUN find $INSTALL_DIR -type f -exec chmod 777 {} \;

#wget --output-document - ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.2-beta-OSX_x86_64.zip | bsdtar xvf - -C /Users/gregpoore/Desktop/AAA_Research_at_UCSD/Oncobiome/Centrifuge --strip-components=1