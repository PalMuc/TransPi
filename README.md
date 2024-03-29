# TransPi - TRanscriptome ANalysiS PIpeline

```text
 _______                                 _____   _
|__   __|                               |  __ \ (_)
   | |     _ __    __ _   _ __    ___   | |__) | _
   | |    |  __|  / _  | |  _ \  / __|  |  ___/ | |
   | |    | |    | (_| | | | | | \__ \  | |     | |
   |_|    |_|     \__,_| |_| |_| |___/  |_|     |_|
 ```

[![Prepint](http://d2538ggaoe6cji.cloudfront.net/sites/default/files/images/favicon.ico)](https://doi.org/10.1101/2021.02.18.431773)[**Preprint**](https://doi.org/10.1101/2021.02.18.431773) &ensp;[![Chat on Gitter](https://img.shields.io/gitter/room/PalMuc/TransPi.svg?colorB=26af64&style=popout)](https://gitter.im/PalMuc/TransPi) &ensp;[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)&ensp;[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)&ensp;[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![release](https://img.shields.io/github/v/release/PalMuc/TransPi?label=release&logo=github)](https://github.com/PalMuc/TransPi/releases/latest)

# Table of contents

* [General info](#General-info)
  * [Pipeline processes](#Pipelie-processes)
  * [Manual](#Manual)
* [Publication](#Publication)
  * [Citation](#Citation)
* [Funding](#Funding)
* [Future work](#Future-work)
* [Issues](#Issues)
  * [Chat](#Chat)

# General info

TransPi – a comprehensive TRanscriptome ANalysiS PIpeline for de novo transcriptome assembly

TransPi is based on the scientific workflow manager [Nextflow](https://www.nextflow.io). It is designed to help researchers get the best reference transcriptome assembly for their organisms of interest. It performs multiple assemblies with different parameters to then get a non-redundant consensus assembly. It also performs other valuable analyses such as quality assessment of the assembly, BUSCO scores, Transdecoder (ORFs), and gene ontologies (Trinotate), etc. All these with minimum input from the user but without losing the potential of a comprehensive analysis.

## Pipeline processes

![TransPi flowchart](https://sync.palmuc.org/index.php/s/nrd3KPnfnz7AipF/preview)

**Figure 1.** TransPi v1.0.0 flowchart showing the various steps and analyses it can performed. For simplicity, this diagram does not show all the connections between the processes. Also, it omits other additional options like the BUSCO distribution and transcriptome filtering with psytrans (see Section 2.6). ORFs=Open reading Frames; HTML=Hypertext Markup Language.

## Manual

TransPi documentation and examples can be found [here](https://palmuc.github.io/TransPi/)

# Publication

Preprint of TransPi including kmer, reads length, and reads quantities tests can be found [here](https://doi.org/10.1101/2021.02.18.431773). Also we tested the pipeline with over 45 samples from different phyla.

TransPi has been peer-reviewed and recommended by Peer Community In Genomics
(https://doi.org/10.24072/pci.genomics.100009)

## Citation

If you use TransPi please cite the peer-reviewed publication:

Rivera-Vicéns, R.E., García-Escudero, CA., Conci, N., Eitel, M., and Wörheide, G. (2021). TransPi – a comprehensive TRanscriptome ANalysiS PIpeline for de novo transcriptome assembly. bioRxiv 2021.02.18.431773; doi: https://doi.org/10.1101/2021.02.18.431773

# Funding

- European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 764840 (ITN IGNITE).

- Advanced Human Capital Program of the National Commission for Scientific and Technological Research (CONICYT)

- Lehre@LMU (project number: W19 F1; Studi forscht@GEO)

- LMU Munich’s Institutional Strategy LMUexcellent within the framework of the German Excellence Initiative

# Future work

- Cloud deployment of the tool

# Issues

We tested TransPi using conda, singularity and docker. However, if you find a problem or get an error please let us know by opening an issue.

## Chat

If you have further questions and need help with TransPi you can chat with us in the [TransPi Gitter chat](https://gitter.im/PalMuc/TransPi)
