FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/bioinf-tools/:$PATH
ENV LANG=C.UTF-8
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

ARG CTE_DIR=/cte
RUN mkdir -p $CTE_DIR/.ci/
COPY .ci/install_dependencies.sh $CTE_DIR/.ci/install_dependencies.sh
RUN $CTE_DIR/.ci/install_dependencies.sh /bioinf-tools

COPY . $CTE_DIR
RUN cd $CTE_DIR \
  && pip3 install tox \
  && cd $CTE_DIR \
  && tox \
  && pip3 install .

CMD cte
