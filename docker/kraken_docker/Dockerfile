# Install Kraken

FROM ubuntu:14.04

ENV INSTALL_DIR /usr/bin/

# Install Kraken
RUN apt-get update
RUN apt-get install -y xz-utils ca-certificates jq wget make g++ sudo

RUN wget --output-document - http://ccb.jhu.edu/software/kraken/dl/kraken-0.10.5-beta.tgz | tar xzf - --directory $INSTALL_DIR  --strip-components=1

ENV PATH ${INSTALL_DIR}:$PATH

RUN cd $INSTALL_DIR && sudo /bin/bash install_kraken.sh $INSTALL_DIR || exit 0

RUN echo 'export PATH="/usr/bin:$PATH"' >> ~/.profile
