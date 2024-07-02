FROM mambaorg/micromamba:1.4.9

LABEL image.author.name "Adolfo Rojas Hidalgo"
LABEL image.author.email "adolfo.rojas@ug.uchile.cl"

USER root
COPY cellranger-8.0.1.tar.gz /opt/cellranger-8.0.1.tar.gz

RUN cd /opt && \
	tar -xzvf cellranger-8.0.1.tar.gz && \
	export PATH=/opt/cellranger-8.0.1:$PATH && \
	ln -s /opt/cellranger-8.0.1/cellranger /usr/bin/cellranger && \
	rm -rf /opt/cellranger-8.0.1.tar.gz # buildkit

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*