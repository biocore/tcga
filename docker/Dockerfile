# Docker version 1.12.0, build 8eab29e
FROM python
COPY . /src
CMD ["python", "/src/cgc_metadata_to_qiime_mapping_file.py"]

# install python requirements
RUN pip install numpy pandas sevenbridges-python
RUN pip install -e git+https://github.com/tqdm/tqdm.git@master#egg=tqdm

# example of command line call:
# docker run gregpoore/cgc2qiime python /src/cgc_metadata_to_qiime_mapping_file.py -f ".fastq" -a AUTH_TOKEN
