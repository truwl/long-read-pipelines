FROM python:3.7.9-slim-buster

MAINTAINER Kiran V Garimella

# install gsutil
RUN apt-get update
RUN apt-get install -y curl git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

RUN git clone https://github.com/broadinstitute/annmas.git
RUN cd annmas && \
    python3 -mvenv venv && \
    . venv/bin/activate && \
    pip install -r dev-requirements.txt && \
    pip install -e .

# activate conda environment
RUN echo "source /annmas/venv/bin/activate" > ~/.bashrc

# We want this to be interactive:
CMD ["/bin/bash"]

