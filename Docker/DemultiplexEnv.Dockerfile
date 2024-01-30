FROM mambaorg/micromamba:1.4.9

LABEL image.author.name "Adolfo Rojas Hidalgo"
LABEL image.author.email "adolfo.rojas@ug.uchile.cl"

COPY --chown=$MAMBA_USER:$MAMBA_USER Demultiplex.yml /tmp/Demultiplex.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER GRCh38_1000G_MAF0.05_ExonFiltered_ChrEncoding.sorted.vcf /tmp/GRCh38_1000G_MAF0.05_ExonFiltered_ChrEncoding.sorted.vcf

RUN micromamba --version

RUN micromamba create -n Demultiplex 

RUN micromamba install -y -n Demultiplex -f /tmp/Demultiplex.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/Demultiplex/bin:$PATH
RUN rm /tmp/Demultiplex.yml
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*