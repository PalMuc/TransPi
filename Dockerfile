FROM continuumio/miniconda3

LABEL authors="Ramon Rivera-Vicens" \
      description="Docker image containing all requirements for TransPi pipeline" \
      version="1.0dev"

RUN apt update; apt install -y gcc bc procps

COPY transpi_env.yml /
RUN conda env create -f /transpi_env.yml  && conda clean -a

ENV PATH /opt/conda/envs/TransPi/bin:$PATH

RUN sed -i 's/base/TransPi/g' ~/.bashrc

RUN wget http://arthropods.eugenes.org/EvidentialGene/other/evigene_older/evigene19may14.tar
RUN tar -xf evigene19may14.tar && rm evigene19may14.tar
ENV PATH /evigene/scripts/prot/:$PATH

RUN mkdir -p /opt/conda/envs/TransPi/lib/python3.6/site-packages/bin && cp /opt/conda/envs/TransPi/bin/skip*.awk /opt/conda/envs/TransPi/lib/python3.6/site-packages/bin/
