FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a

# install flair
RUN wget https://github.com/BrooksLabUCSC/flair/archive/v1.4.tar.gz \
        && tar xf v1.4.tar.gz

# set path
ENV PATH=/opt/conda/envs/lr-flair/bin/:/flair-1.4/:${PATH}
