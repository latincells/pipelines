FROM mambaorg/micromamba:1.4.9

LABEL image.author.name "Adolfo Rojas Hidalgo"
LABEL image.author.email "adolfo.rojas@ug.uchile.cl"

COPY --chown=$MAMBA_USER:$MAMBA_USER LatinCells.yml /tmp/LatinCells.yml

RUN micromamba --version

RUN micromamba create -n SingleCell 

RUN micromamba install -y -n SingleCell -f /tmp/LatinCells.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/SingleCell/bin:$PATH
RUN R -e 'install.packages("Matrix", type = "source", repos="https://cran.dcc.uchile.cl")'
RUN R -e 'install.packages("irlba", type = "source", repos="https://cran.dcc.uchile.cl")'
#RUN R -e 'install.packages("devtools", type = "source", repos="https://cran.dcc.uchile.cl")'
#RUN R -e 'install.packages("NMF",repos="https://cran.dcc.uchile.cl")'
RUN R -e 'remotes::install_github("satijalab/seurat-wrappers", upgrade=F)'
RUN R -e 'remotes::install_github("saezlab/liana", upgrade=F)'
RUN R -e 'devtools::install_github("immunogenomics/harmony",upgrade ="never")'
RUN pip install umap-learn==0.5.3
RUN R -e 'devtools::install_github("immunogenomics/presto",upgrade ="never")'
RUN R -e 'devtools::install_github("jinworks/CellChat",upgrade ="never")'
#RUN R -e 'BiocManager::install("clusterExperiment",update = FALSE)'
RUN rm /tmp/LatinCells.yml
USER root
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*