FROM continuumio/miniconda3

LABEL authors="Ramon Rivera-Vicens" \
      description="Docker image containing all requirements for TransPi pipeline"

RUN apt update; apt install -y gcc bc procps

COPY transpi_env.yml /
RUN conda env create -f /transpi_env.yml  && conda clean -a

ENV PATH /opt/conda/envs/TransPi/bin:$PATH
