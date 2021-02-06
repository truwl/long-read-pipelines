FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# install gsutil
RUN apt update
RUN apt install -y curl git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash

RUN git clone https://github.com/broadinstitute/annmas.git
RUN cd annmas && \
    python3 -mvenv venv && \
    . venv/bin/activate && \
    pip install -r dev-requirements.txt && \
    pip install -e .

# activate conda environment
RUN echo "source /annmas/venv/bin/activate" > ~/.bashrc
