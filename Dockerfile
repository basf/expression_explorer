FROM rocker/shiny:latest

RUN apt-get update -qq && apt-get install -y --no-install-recommends libxml2-dev

RUN  R -e "install.packages(c('ggplot2','ggrepel','gridExtra','reshape2','BiocManager'))"
RUN  R -e "install.packages(c('dplyr','ape','DT'))"
RUN  R -e "install.packages('devtools')"
RUN  R -e "library('devtools'); devtools::install_github('CosteaPaul/r-cytoscape.js')"
RUN  R -e "install.packages(c('shinyBS','VennDiagram','shinyWidgets'))"
RUN  R -e "remotes::install_github('Marlin-Na/trewjb')"
RUN  R -e "BiocManager::install('DESeq2')"
RUN  R -e "install.packages('readr')"

ADD app /srv/shiny-server/expression_explorer
