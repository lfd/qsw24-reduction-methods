FROM ubuntu:22.04

LABEL maintainer="lfd"

WORKDIR /repro

COPY ./code/ ./

RUN apt update \
    && apt install -y build-essential \
    && apt install -y wget \
    && apt clean \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
	
ENV PATH=$CONDA_DIR/bin:$PATH

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" >> /etc/apt/sources.list 
RUN apt-get update
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get install -y r-base

RUN R -e "install.packages('ggplot2', dependencies=TRUE)"
RUN R -e "install.packages('latex2exp', dependencies=TRUE)"
RUN R -e "install.packages('ggsci', dependencies=TRUE)"
RUN R -e "install.packages('gridExtra', dependencies=TRUE)"
RUN R -e "install.packages('grid', dependencies=TRUE)"
RUN R -e "install.packages('patchwork', dependencies=TRUE)"
RUN R -e "install.packages('tikzDevice', dependencies=TRUE)"
RUN R -e "install.packages('stringr', dependencies=TRUE)"
RUN R -e "install.packages('scales', dependencies=TRUE)"
RUN R -e "install.packages('dplyr', dependencies=TRUE)"
RUN R -e "install.packages('tidyr', dependencies=TRUE)"
RUN R -e "install.packages('readr', dependencies=TRUE)"

RUN apt install -y texlive-full

RUN conda create -n quark_install python=3.10

SHELL ["conda", "run", "-n", "quark_install", "/bin/bash", "-c"]

RUN conda config --add channels conda-forge
RUN conda update -n base -c defaults conda
RUN conda install -c dlr-sc quark=0.1
RUN conda install -c conda-forge pathos
RUN conda install pandas
RUN conda install scipy
RUN pip install pyqubo
RUN pip install qiskit==0.41.1

CMD ["/bin/bash"]
