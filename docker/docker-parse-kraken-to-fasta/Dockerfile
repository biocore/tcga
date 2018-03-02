# Creates Dockerfile for kraken_to_fasta.py
FROM python:3
MAINTAINER Evguenia Kopylova <jenya.kopylov@gmail.com>

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

RUN echo nameserver 8.8.8.8 > /etc/resolv.conf && pip install numpy
RUN echo nameserver 8.8.8.8 > /etc/resolv.conf && pip install click
RUN echo nameserver 8.8.8.8 > /etc/resolv.conf && pip install scikit-bio
COPY kraken_to_fasta.py /usr/src/app
COPY test_kraken_to_fasta.py /usr/src/app

ADD kraken_to_fasta.py /usr/src/app/kraken_to_fasta.py

CMD ["python", "/usr/src/app/kraken_to_fasta.py"]
