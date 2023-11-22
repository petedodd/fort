# docker build -t tbss .
# docker tag tbss swtoolscr.azurecr.io/tbstatisticalserver:beta
# docker push swtoolscr.azurecr.io/tbstatisticalserver:beta
# docker run -d -p 8080:8080 --name tbss swtoolscr.azurecr.io/tbstatisticalserver:beta

FROM rstudio/plumber

RUN Rscript -e "install.packages('devtools')"

ARG GIT_BRANCH=main

WORKDIR  /R

RUN git clone -b $GIT_BRANCH --single-branch https://github.com/petedodd/fort.git fort
RUN Rscript -e "devtools::install('fort')"

COPY startup.R .

EXPOSE 8000

CMD ["Rscript startup.R"]