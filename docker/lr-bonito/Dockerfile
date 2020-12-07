FROM nvcr.io/nvidia/pytorch:20.11-py3

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y wget gcc make git
RUN pip install --use-feature=2020-resolver ont-bonito==0.3.1
