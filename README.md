# fort <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/petedodd/fort/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/petedodd/fort/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

R package to forecast TB notifications, incidence, mortality, and prevalence

## License

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Usage

While not public, the easiest approach is to clone the repository, and then use:

```R
devtools::install('fort')
```

assuming the working directory is one level up from the repository clone.

## API

A plumber API is now defined in inst/plumber.

This can be run using:

```R
library(plumber)
root <- plumb_api(package='fort',name='tbstatisticalserver')
pr_run(root)
```

The server can be tested using the json data in the same folder:

```bash
curl -X "POST" "http://127.0.0.1:8021/projection" -H 'Content-Type: application/json' -d @projection.json
```

## Docker

fort can be packaged as a Docker image. The relevant files are in the `/docker` folder.

The files are cloned (indirectly) from the GitHub repo, for historial reasons.

To build the Docker image, use:

```bash
docker build -t fort .
```

To (optionally) push to the Docker image to the Avenir Health container registry (permissions required), use:

```bash
docker tag fort swtoolscr.azurecr.io/tbstatisticalserver:beta
docker push swtoolscr.azurecr.io/tbstatisticalserver:beta
```

To launch a Docker container from the image, use:

```bash
docker run -d -p 8080:8080 --name fort fort
```
