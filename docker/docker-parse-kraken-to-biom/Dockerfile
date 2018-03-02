# Creates Dockerfile for parse_kraken_to_biom.py
FROM python:3
MAINTAINER Evguenia Kopylova <jenya.kopylov@gmail.com>

RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app

COPY requirements.txt /usr/src/app/
RUN pip install --no-cache-dir numpy h5py
RUN pip install --no-cache-dir -r requirements.txt
COPY parse_kraken_to_biom.py /usr/src/app
COPY test_parse_kraken_to_biom.py /usr/src/app

ADD parse_kraken_to_biom.py /usr/src/app/parse_kraken_to_biom.py

CMD ["python", "/usr/src/app/parse_kraken_to_biom.py”]
