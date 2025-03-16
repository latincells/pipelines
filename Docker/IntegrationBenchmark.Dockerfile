FROM mambaorg/micromamba:1.4.9

LABEL image.author.name "Adolfo Rojas Hidalgo"
LABEL image.author.email "adolfo.rojas@ug.uchile.cl"

COPY --chown=$MAMBA_USER:$MAMBA_USER LatinCells.IntengrationBenchmark.yml /tmp/LatinCells.IntengrationBenchmark.yml

RUN micromamba --version

RUN micromamba create -n LatinCells 

RUN micromamba install -y -n LatinCells -f /tmp/LatinCells.IntengrationBenchmark.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/LatinCells/bin:$PATH
RUN pip install scib-metrics
RUN pip install datashader
RUN pip install colorcet
RUN pip install bokeh
RUN pip install holoviews
RUN pip install scanorama
RUN pip install harmony-pytorch
RUN pip install scvi-tools
RUN rm /tmp/LatinCells.IntengrationBenchmark.yml
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*