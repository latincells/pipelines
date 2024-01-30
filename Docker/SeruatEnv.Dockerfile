FROM mambaorg/micromamba:1.4.9

LABEL image.author.name "Adolfo Rojas Hidalgo"
LABEL image.author.email "adolfo.rojas@ug.uchile.cl"

COPY --chown=$MAMBA_USER:$MAMBA_USER LatinCells.yml /tmp/LatinCells.yml

RUN micromamba --version

RUN micromamba create -n SingleCell 

RUN micromamba install -y -n SingleCell -f /tmp/LatinCells.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/SingleCell/bin:$PATH
RUN rm /tmp/LatinCells.yml
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*