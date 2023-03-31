FROM debian:10

ENV DEBIAN_FRONTEND=noninteractive

# Install system libraries required for python and R installations
RUN apt-get update --fix-missing && apt-get install -y --no-install-recommends build-essential apt-utils ca-certificates zlib1g-dev gfortran locales libxml2-dev libopenmpi-dev libcurl4-openssl-dev libssl-dev libzmq3-dev libreadline6-dev xorg-dev libcairo-dev libpango1.0-dev libbz2-dev liblzma-dev libffi-dev libsqlite3-dev libhdf5-dev libjpeg-dev libblas-dev liblapack-dev libpcre2-dev libgit2-dev libgmp-dev
RUN apt-get install -y --no-install-recommends libpcre2-dev libigraph0-dev
# Install common linux tools
RUN apt-get update && apt-get install -y --no-install-recommends wget curl htop less nano vim emacs git

# Configure default locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8



# Download and compile R from source
WORKDIR /app/R
RUN wget https://cran.rstudio.com/src/base/R-4/R-4.0.2.tar.gz
RUN tar xvfz R-4.0.2.tar.gz && rm R-4.0.2.tar.gz

WORKDIR /app/R/R-4.0.2
RUN ./configure --enable-R-shlib --with-cairo --with-libpng --prefix=/app/R/
RUN make && make install
ENV PATH="/app/R/bin:${PATH}"
ENV LD_LIBRARY_PATH="/app/R/lib/R/lib:${LD_LIBRARY_PATH}"

RUN Rscript -e "update.packages(ask=FALSE, repos='https://ftp.gwdg.de/pub/misc/cran/')"
RUN Rscript -e "install.packages(c('devtools', 'gam', 'RColorBrewer', 'UpSetR', 'sctransform', 'BiocManager'), repos='https://ftp.gwdg.de/pub/misc/cran/')"
RUN Rscript -e "devtools::install_github(c('patzaw/neo2R', 'patzaw/BED'))"
RUN Rscript -e "BiocManager::install(c('scran', 'dorothea', 'MAST', 'limma', 'viper', 'Seurat','ComplexHeatmap','slingshot','clusterExperiment'), version = '3.12')"


#RUN Rscript -e 'writeLines(capture.output(sessionInfo()), "../../package_versions_r.txt")' --default-packages=scran,RColorBrewer,slingshot,monocle,gam,clusterExperiment,ggplot2,plyr



# Download and compile python from source
WORKDIR /app/python3
RUN wget https://www.python.org/ftp/python/3.8.7/Python-3.8.7.tgz
RUN tar zxfv Python-3.8.7.tgz && rm Python-3.8.7.tgz

WORKDIR /app/python3/Python-3.8.7
RUN ./configure --with-lto --prefix=/app/python3/
RUN make && make install

WORKDIR /app/python3
RUN rm -rf /opt/python3/Python-3.8.7
RUN ln -s /app/python3/bin/python3 /app/python3/bin/python
RUN ln -s /app/python3/bin/pip3 /app/python3/bin/pip
ENV PATH="/app/python3/bin:${PATH}"

COPY python-packages.txt /app/python3/python-packages.txt
RUN pip install --no-cache-dir -U pip wheel setuptools cmake
RUN pip install --no-cache-dir -r /app/python3/python-packages.txt
RUN jupyter contrib nbextension install --system
RUN jupyter nbextension enable --py widgetsnbextension

RUN pip freeze > ../../package_versions_py.txt

COPY .bashrc_docker /root/.bashrc
COPY .profile_docker /root/.profile

RUN apt-get clean -y && apt-get autoremove -y
