# Use a Rocker image as the base image, specifying the version of R you need
FROM rocker/r-ver:4.3.2

# Install Linux dependencies (you may need more depending on your R packages)
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    pandoc \
    pandoc-citeproc \
    libuv1-dev \
    zlib1g-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev


# Install the renv package globally
RUN R -e "install.packages('renv')"

# Set the working directory to the container
WORKDIR /root

# Copy entire project directory into the container
COPY . /root

# Change default location of cache to the container
# (cache means that renv will not be restored from scratch every time)
RUN mkdir renv/.cache
ENV RENV_PATHS_CACHE renv/.cache

# Restore the project library using renv
RUN R -e "renv::restore()"

# Render Rmarkdown file
# (using CMD instead of RUN renders this at run time, not build time)
# (using RUN would not work as environment not set)
CMD ["R", "-e", "rmarkdown::render('code/final_analytics/paper-from-code.Rmd')"]
