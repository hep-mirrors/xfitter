FROM gitlab-registry.cern.ch/fitters/xfitter:with_docker
WORKDIR /home/xfitter
COPY . .
RUN ./tools/install-xfitter deps
RUN . ./setup.sh
RUN ./make.sh install
CMD ["cd" , "/home/xfitter"]
