FROM quay.io/bitnami/python:3.7.11

RUN mkdir /gex_analysis

WORKDIR /gex_analysis

COPY . .

# Install GEX_ANALYSIS and its dependencies
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r requirements.txt
