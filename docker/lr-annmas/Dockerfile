FROM python:3.7.9-slim-buster

MAINTAINER Kiran V Garimella

# install gsutil
RUN apt-get update
RUN apt-get install -y curl git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash

RUN git clone https://github.com/broadinstitute/annmas.git
RUN cd annmas && git checkout jts_sirv_alevin_input  
RUN cd annmas && \
    python3 -mvenv venv && \
    . venv/bin/activate && \
    pip install -r dev-requirements.txt && \
    pip install -e .

# activate conda environment
RUN echo "source /annmas/venv/bin/activate" > ~/.bashrc

# We want this to be interactive:
CMD ["/bin/bash"]

