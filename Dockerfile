#Pick you favourite shiny image
FROM shiny_image:latest

RUN apt-get update -qq && apt-get install -y --no-install-recommends libtiff-dev openssh-client

RUN  R -e "remotes::install_cran(c('shiny', 'rmarkdown','ggplot2','XLConnect','ggrepel','gridExtra','imager','reshape2','BiocManager','plotly'), repos='https://cran.rstudio.com/')" && \
      R -e "BiocManager::install('DESeq2')"&& \
      cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
      rm -rf /var/lib/apt/lists/*

RUN  R -e "install.packages(c('dplyr','ape','DT'), repos='https://cran.rstudio.com/')"&& \
      cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/

RUN  R -e "install.packages('devtools', repos='https://cran.rstudio.com/')"
RUN  R -e "library('devtools'); devtools::install_github('CosteaPaul/r-cytoscape.js')"
RUN  R -e "install.packages(c('shinyBS','VennDiagram','shinyWidgets'), repos='https://cran.rstudio.com/')"
RUN  R -e "remotes::install_github('Marlin-Na/trewjb')"
# If you need database-drivers, use shiny_dbdrivers:latest

ARG PROJECT_NAME

## to add additional system libs or packages uncomment the next lines

## Install additional system libs (ubuntu!)
# ARG DEBIAN_FRONTEND=noninteractive

## Install additional R packages
# RUN Rscript -e 'remotes::install_cran(c("<Package1>", "<Package2>"), upgrade = "always", repos = "https://cran.rstudio.com/")'

ADD app /srv/shiny-server/$PROJECT_NAME
USER 999:999
