FROM quay.io/jitesoft/ubuntu:18.04

RUN mkdir /gex_analysis

WORKDIR /gex_analysis

COPY . .

# Installation - Dep
RUN apt-get update && apt-get install -y python3-setuptools libssl-dev libffi-dev build-essential

# Install GEX_ANALYSIS and its dependencies
RUN pip3 install -r requirements.txt
