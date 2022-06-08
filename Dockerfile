FROM bgruening/galaxy-stable:20.09
USER galaxy
COPY geec_galaxy /galaxy-central/tools/geec_galaxy
RUN rm /galaxy-central/config/tool_conf.xml.main
RUN rm /galaxy-central/config/tool_conf.xml.sample
COPY tool_conf.xml /galaxy-central/config/tool_conf.xml
#COPY galaxy.ini /galaxy-central/config/galaxy.ini
USER root