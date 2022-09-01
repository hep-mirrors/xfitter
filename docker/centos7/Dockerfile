FROM rootproject/root:6.26.06-centos7
LABEL maintainer.name="xFitter Developers team"
LABEL maintainer.email="xfitter-users@googlegroups.com"

USER root
ENV LANG=C.UTF-8
ENV XFITTER_VERSION master
ENV XFITTER_INSTALL_DIR /opt/xfitter
WORKDIR $XFITTER_INSTALL_DIR

RUN yum -y install yum-plugin-ovl which yaml-cpp-devel libyaml-devel gsl-devel cmake automake libtool wget cmake \
    && rm -rf /packages /var/lib/apt/lists/*

RUN wget https://gitlab.cern.ch/fitters/xfitter/-/raw/docker/tools/install-xfitter \
    && chmod +x install-xfitter \
    && ./install-xfitter ${XFITTER_VERSION} \
    && chmod -R 755 ${XFITTER_INSTALL_DIR} \
    && rm -rf install-xfitter install.log run deps/lhapdf/share/LHAPDF/NNPDF30_nlo_as_0118

COPY docker-entrypoint.sh /
RUN chmod +x /docker-entrypoint.sh

VOLUME [ "/run" ]
VOLUME [ "/data" ]
VOLUME [ "/pdfdata" ]

WORKDIR /run
ENTRYPOINT [ "/docker-entrypoint.sh" ]

CMD [ "xfitter", "&&", "xfitter-draw", "output" ]