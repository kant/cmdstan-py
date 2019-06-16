# cmdstan-py
This repo also includes `cmdstan.py`, a thin wrapper around cmdstan, and BayesGLM, a command-line tool to run Bayessian Regression on an expression dataset

#### Requires:
```
numpy, pandas, Cef_tools, loompy
```
## To install cmdstan do

```
wget https://github.com/stan-dev/cmdstan/releases/download/v2.12.0/cmdstan-2.12.0.tar.gz
tar -xzf cmdstan-2.12.0.tar.gz
cd cmdstan-2.12.0
make build -j8
```

## Usage Example
```
./BayesGLM.py -i mycef_filepath.cef -s /path/to/cmdstan-2.12.0/ -p 30
```
or better, use `nohup` since it will run for long
```
nohup ./BayesGLM.py -i mycef_filepath.cef -s /path/to/cmdstan-2.12.0/ -p 30 &
```
