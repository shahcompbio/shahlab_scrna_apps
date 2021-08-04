FROM quay.io/ubuntu:18.04

RUN mkdir /gex_analysis

WORKDIR /gex_analysis

# Installation - Dep
RUN apt-get update && apt-get install -y python3-setuptools libssl-dev libffi-dev build-essential

# Install GEX_ANALYSIS and its dependencies
COPY ../../apps/gex_analysis .
RUN pip3 install -r requirements.txt
