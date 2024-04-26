FROM mambaorg/micromamba:1.4.9

LABEL image.author.name "Adolfo Rojas Hidalgo"
LABEL image.author.email "adolfo.rojas@ug.uchile.cl"

COPY --chown=$MAMBA_USER:$MAMBA_USER Celltypist.yml /tmp/Celltypist.yml

RUN micromamba --version

RUN micromamba create -n Celltypist 

RUN micromamba install -y -n Celltypist -f /tmp/Celltypist.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/Celltypist/bin:$PATH
RUN rm /tmp/Celltypist.yml
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*