# Trimmed Harrell-Davis quantile estimator based on the highest density interval of the given width

The repository contains auxiliary data for the paper "Trimmed Harrell-Davis quantile estimator based on the highest density interval of the given width":

* The source code of the paper: `thdqe.Rmd`, `preamble.tex`, `before_body.tex`, `references.bib`, `*.R`
* Calculated data for Simulation 2: `efficiency.csv`
* The reference implementation of the suggested estimator: `thdqe.R`

The paper was published by Taylor & Francis in Communications in Statistics — Simulation and Computation on 17 March 2022, available online: https://www.tandfonline.com/10.1080/03610918.2022.2050396

You can cite it as follows:

> Andrey Akinshin (2022) Trimmed Harrell-Davis quantile estimator based on the highest density interval of the given width, Communications in Statistics - Simulation and Computation, DOI: [10.1080/03610918.2022.2050396](https://doi.org/10.1080/03610918.2022.2050396)

BiBTeX reference:

```bib
@article{akinshin2022thdqe,
  author = {Akinshin, Andrey},
  title = {Trimmed Harrell-Davis quantile estimator based on the highest density interval of the given width},
  journal = {Communications in Statistics - Simulation and Computation},
  pages = {1-11},
  year = {2022},
  month = {3},
  publisher = {Taylor & Francis},
  doi = {10.1080/03610918.2022.2050396},
  URL = {https://www.tandfonline.com/doi/abs/10.1080/03610918.2022.2050396},
  eprint = {https://www.tandfonline.com/doi/pdf/10.1080/03610918.2022.2050396},
  abstract = {Traditional quantile estimators that are based on one or two order statistics are a common way to estimate distribution quantiles based on the given samples. These estimators are robust, but their statistical efficiency is not always good enough. A more efficient alternative is the Harrell-Davis quantile estimator which uses a weighted sum of all order statistics. Whereas this approach provides more accurate estimations for the light-tailed distributions, it’s not robust. To be able to customize the tradeoff between statistical efficiency and robustness, we could consider a trimmed modification of the Harrell-Davis quantile estimator. In this approach, we discard order statistics with low weights according to the highest density interval of the beta distribution.}
}
```

arXiv preprint: [arXiv:2111.11776](https://arxiv.org/abs/2111.11776) **[stat.ME]**