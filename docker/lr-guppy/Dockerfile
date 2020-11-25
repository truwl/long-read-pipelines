FROM nvidia/cuda:10.2-base-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y wget lsb-release
RUN wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add -
RUN echo "deb http://mirror.oxfordnanoportal.com/apt xenial-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list

RUN apt-get update
RUN apt-get install -yq ont-guppy=4.2.2-1~xenial -V --no-install-recommends
RUN apt-get install -yq git gcc make python3-pip virtualenv

ENV PATH=/taiyaki/bin/:${PATH}

RUN git clone https://github.com/nanoporetech/taiyaki && \
        cd taiyaki && \
        make install && \
        . venv/bin/activate

RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && \
        tar xjf minimap2-2.17_x64-linux.tar.bz2 && \
        rm minimap2-2.17_x64-linux.tar.bz2 && \
        cd minimap2-2.17_x64-linux && \
        cp minimap2 /usr/local/bin/
