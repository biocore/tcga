# Creates Dockerfile for parse_kraken_to_biom.py
FROM python:3
MAINTAINER Evguenia Kopylova <jenya.kopylov@gmail.com>

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

COPY requirements.txt /usr/src/app/
RUN pip install --no-cache-dir numpy
RUN pip install --no-cache-dir -r requirements.txt
COPY . /usr/src/app

ADD parse_kraken_to_biom.py /usr/src/app/parse_kraken_to_biom.py
