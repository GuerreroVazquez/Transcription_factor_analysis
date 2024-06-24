FROM jupyter/datascience-notebook:x86_64-ubuntu-22.04

USER root
RUN apt update
RUN pip install anndata igraph scanpy decoupler liana[full] adjustText pybiomart wordcloud
RUN apt install -y graphviz
RUN pip install seaborn --upgrade

USER amore
WORKDIR /home/amore
