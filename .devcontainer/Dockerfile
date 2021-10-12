# Get official docker image fo julia 
FROM julia

# set username as `root`
USER root

# update package list
RUN apt-get update

# install git and x11-apps for x11 test
RUN apt-get install -y git x11-apps

# Required package for GR backend of Plots
# https://gr-framework.org/julia.html#installation
RUN apt-get install -y libxt6 libxrender1 libxext6 libgl1-mesa-glx libqt5widgets5 openssh-client socat

RUN apt-get -yq install git
RUN git config --global user.name "KoheiToyoda"
RUN git config --global user.email "toyoda@chiba-u.jp"

# Install `Plots.jl and precompile it`
RUN julia -e 'using Pkg; Pkg.add(["Plots","BenchmarkTools"])'
RUN julia -e 'ENV["CONDA_JL_VERSION"]=3'
RUN julia -e 'using Pkg; Pkg.add(["LinearAlgebra","Parameters","Polynomials"])'
RUN julia -e 'using Pkg; Pkg.add(["CSV","DataFrames"])'
RUN julia -e 'using Pkg; Pkg.add(["ColorSchemes","PyCall","Pluto"])'
RUN julia -e 'using Pkg; Pkg.add(url = "https://github.com/KoheiToyoda/KTOptical")'
# RUN julia -e 'using Pkg; Pkg.add(["Plotly","PlotlyJS"])'
RUN julia -e 'using Pkg; Pkg.add(["PyPlot"])'
RUN julia -e 'using Pkg; Pkg.add(["JLD2"])'
RUN julia -e 'using Pkg; Pkg.precompile()'
RUN julia -e 'using Plots;'

# Configuration for x11 forwarding
ENV DISPLAY host.docker.internal:0.0
ENV XDG_RUNTIME_DIR /workspaces/Julia/.tmp
ENV LIBGL_ALWAYS_INDIRECT 1


EXPOSE 1234
VOLUME /julia

ENV JULIA_NUM_THREADS = 100
WORKDIR /julia

ENTRYPOINT [ "julia" ]