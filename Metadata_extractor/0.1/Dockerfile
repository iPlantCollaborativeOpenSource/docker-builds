From ubuntu:14.04.3
MAINTAINER HUANIAN ZHANG <fantasyzhn@email.arizona.edu>
LABEL "This Dockerfile is for Metadata_extractor:ME_v01"

RUN apt-get update && apt-get install -y git make gcc wget python-dev python2.7 python-numpy python-pip matplotlib
RUN pip install --no-deps astropy
ENV BINPATH /usr/bin
ADD extract_metadata.py $BINPATH
RUN chmod +x $BINPATH/extract_metadata.py

ENTRYPOINT ["extract_metadata.py"]
