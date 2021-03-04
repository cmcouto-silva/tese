FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y locales
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && locale-gen

ENV LANG C.UTF-8
ENV LANGUAGE C.UTF-8
ENV LC_CTYPE "C.UTF-8"
ENV LC_NUMERIC "C.UTF-8"
ENV LC_TIME "C.UTF-8"
ENV LC_COLLATE "C.UTF-8"
ENV LC_MONETARY "C.UTF-8"
ENV LC_MESSAGES "C.UTF-8"
ENV LC_PAPER "C.UTF-8"
ENV LC_NAME "C.UTF-8"
ENV LC_ADDRESS "C.UTF-8"
ENV LC_TELEPHONE "C.UTF-8"
ENV LC_MEASUREMENT "C.UTF-8"
ENV LC_IDENTIFICATION "C.UTF-8"
ENV LC_ALL C.UTF-8

WORKDIR /app

ADD config config/

RUN bash config/install.sh
RUN bash config/requirements.sh
RUN Rscript config/requirements.R
RUN pip3 install -r config/requirements.txt

CMD ["/bin/bash"]
