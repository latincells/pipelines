FROM mambaorg/micromamba:1.4.9

LABEL image.author.name "Adolfo Rojas Hidalgo"
LABEL image.author.email "adolfo.rojas@ug.uchile.cl"

COPY --chown=$MAMBA_USER:$MAMBA_USER LatinCells.Hackaton.yml /tmp/LatinCells.Hackaton.yml

RUN micromamba --version

RUN micromamba create -n LatinCells 

RUN micromamba install -y -n LatinCells -f /tmp/LatinCells.Hackaton.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/LatinCells/bin:$PATH
RUN R -e 'install.packages("Matrix", type = "source", repos="https://cran.dcc.uchile.cl")'
RUN R -e 'install.packages("irlba", type = "source", repos="https://cran.dcc.uchile.cl")'
RUN R -e 'remotes::install_github("satijalab/seurat-wrappers", upgrade=F)'
RUN R -e 'remotes::install_github("saezlab/liana", upgrade=F)'
RUN R -e 'devtools::install_github("immunogenomics/harmony",upgrade ="never")'
RUN pip install umap-learn==0.5.3
RUN R -e 'install.packages("NMF",repos="https://cran.dcc.uchile.cl")'
RUN R -e 'devtools::install_github("jinworks/CellChat",upgrade ="never")'
RUN R -e 'devtools::install_github("smorabit/hdWGCNA", ref="dev",upgrade ="never")'
#RUN R -e 'devtools::install_github("satijalab/seurat-data")'
RUN rm /tmp/LatinCells.Hackaton.yml
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*