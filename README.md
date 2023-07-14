# fort <img src="man/figures/logo.png" align="right" height="139" />

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


## TODO

- fix mortality offset
- introduce hyperparameter for IP model
