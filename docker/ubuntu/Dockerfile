FROM rootproject/root:6.24.06-ubuntu20.04
LABEL maintainer.name="xFitter Developers team"
LABEL maintainer.email="xfitter-users@googlegroups.com"

USER root
ENV XFITTER_VERSION master
ENV XFITTER_INSTALL_DIR /opt/xfitter
WORKDIR $XFITTER_INSTALL_DIR

RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends wget cmake automake libtool libyaml-dev libyaml-cpp-dev \
    && rm -rf /packages /var/lib/apt/lists/*

RUN wget https://gitlab.cern.ch/fitters/xfitter/-/raw/master/tools/install-xfitter \
    && chmod +x install-xfitter \
    && ./install-xfitter ${XFITTER_VERSION} \
    && chmod -R 755 ${XFITTER_INSTALL_DIR} \
    && rm -rf install-xfitter install.log run deps/lhapdf/share/LHAPDF/NNPDF30_nlo_as_0118

#COPY docker-entrypoint.sh /
#RUN chmod +x /docker-entrypoint.sh

#VOLUME [ "/run" ]
#VOLUME [ "/data" ]
#VOLUME [ "/pdfdata" ]

WORKDIR /run
#ENTRYPOINT [ "/docker-entrypoint.sh" ]
CMD [ "xfitter", "&&", "xfitter-draw", "output" ]