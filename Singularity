Bootstrap:docker
From:continuumio/miniconda3

%labels
    MAINTAINER Ramon Rivera-Vicens
    DESCRIPTION Singularity image containing all requirements for TransPi pipeline
    VERSION 1.0dev
    
%environment
    PATH=/opt/conda/envs/TransPi/bin:$PATH
    export PATH

%files
    transpi_env.yml /

%post
    apt update; apt install -y gcc bc procps
    /opt/conda/bin/conda env create -f /transpi_env.yml
    /opt/conda/bin/conda clean -a
