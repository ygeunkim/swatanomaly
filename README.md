
# `swatanomaly`

Various anomaly detection algorithms studied in *NSR anomaly detection
project*

  - NND (baseline-method): Yun et al. (2018)
  - K-L divergence: Cho et al. (2019)
  - SOM

## Installation

``` r
devtools::install_github("ygeunkim/swatanomaly")
```

``` r
library(swatanomaly)
```

## Nearest-Neighbor Distance Algorithm

This algorithm follows k-fold cross-validation construction.

### Windowed NNS

| 123456789 | 123456789 | 123456789 | 123456789 | 123456789 |
| :-------: | :-------: | :-------: | :-------: | :-------: |
| ooooooooo | ooooooooo | ooooooooo | ooooooooo | ooooooooo |
| ooooooooo | ooooooooo | ooooooooo | ooooooooo | ooooooooo |
| ooooooooo | ooooooooo | ooooooooo | ooooooooo | ooooooooo |

As k-fold CV, one fold is validation block and the remainings are
training. NND is computed versus training. Consider the first block.
Then we should compute NND of 9 observations, respectively.

For each ![i
\\leftarrow 1](https://latex.codecogs.com/png.latex?i%20%5Cleftarrow%201
"i \\leftarrow 1") to w(size of the block),

1.  Compute **Euclidean distance** versus every point in the training
    block.
2.  Find minimum value of the distance value.
3.  It is the NND of this observation.

Repeat this process for every block. Here, Euclidean distance for
multivariate series is computed by

  
![d(\\mathbf{x}, \\mathbf{y}) = \\sqrt{\\sum\_j^p (x\_j -
y\_j)^2}](https://latex.codecogs.com/png.latex?d%28%5Cmathbf%7Bx%7D%2C%20%5Cmathbf%7By%7D%29%20%3D%20%5Csqrt%7B%5Csum_j%5Ep%20%28x_j%20-%20y_j%29%5E2%7D
"d(\\mathbf{x}, \\mathbf{y}) = \\sqrt{\\sum_j^p (x_j - y_j)^2}")  

See `euc_nnd()` and `nns_cpp()`.

### Setting threshold

1.  arbitrary threshold
2.  right tail value of NND pdf

`euc_pdf()` produces this NND pdf. The procedure is similar to the above
windowed NNS with smaller partition, but different input.

  - `nns_cpp()` is provided window size (`win`), while
  - `euc_pdf()` partition number (`partition`)

Using this pdf, we can compute its quantile or all at once using
`nnd_thr()`.

### Anomaly detection

Given the result of windowed NNS, compare with the threshold. Use
`detect_nnd()`. It returns `TRUE` (anomaly) and `FALSE` (normal).
Sometimes we want the other labels. In this case, `detect_nndvec` can
give any label you want (`label`). For example, `-1` and `1`.

## Kullback-Leibler Divergence Algorithm

This algorithm implements neighboring-window method. Using
`stats::density()`, we first estimate each probability mass. Then we can
compute K-L divergence from ![p](https://latex.codecogs.com/png.latex?p
"p") to ![q](https://latex.codecogs.com/png.latex?q "q") by

  
![D\_{KL} = \\sum\_{x \\in \\mathcal{X}} p(x) \\ln
\\frac{p(x)}{q(x)}](https://latex.codecogs.com/png.latex?D_%7BKL%7D%20%3D%20%5Csum_%7Bx%20%5Cin%20%5Cmathcal%7BX%7D%7D%20p%28x%29%20%5Cln%20%5Cfrac%7Bp%28x%29%7D%7Bq%28x%29%7D
"D_{KL} = \\sum_{x \\in \\mathcal{X}} p(x) \\ln \\frac{p(x)}{q(x)}")  

where
![\\mathcal{X}](https://latex.codecogs.com/png.latex?%5Cmathcal%7BX%7D
"\\mathcal{X}") is the support of
![p](https://latex.codecogs.com/png.latex?p "p")
(![f\_1](https://latex.codecogs.com/png.latex?f_1 "f_1")).

The algorithm uses a threshold
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\\lambda")
and check if ![D \<
\\lambda](https://latex.codecogs.com/png.latex?D%20%3C%20%5Clambda
"D \< \\lambda"). If this holds, two windows can be said that they are
derived from the same gaussian distribution.

### Fixed lambda algorithm

  - Data: univariate series of size
    ![n](https://latex.codecogs.com/png.latex?n "n")
  - Input: window size ![win](https://latex.codecogs.com/png.latex?win
    "win"), jump size for sliding window
    ![jump](https://latex.codecogs.com/png.latex?jump "jump"), threshold
    ![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda
    "\\lambda")

Note that the number of window is

  
![\\frac{n - win}{jump} + 1 \\equiv
w](https://latex.codecogs.com/png.latex?%5Cfrac%7Bn%20-%20win%7D%7Bjump%7D%20%2B%201%20%5Cequiv%20w
"\\frac{n - win}{jump} + 1 \\equiv w")  

1.  For ![i
    \\leftarrow 1](https://latex.codecogs.com/png.latex?i%20%5Cleftarrow%201
    "i \\leftarrow 1") to (w - 1) do
    1.  if ![D \<
        \\lambda](https://latex.codecogs.com/png.latex?D%20%3C%20%5Clambda
        "D \< \\lambda") then
        1.  ![f\_1](https://latex.codecogs.com/png.latex?f_1 "f_1") be
            the pdf of ![i](https://latex.codecogs.com/png.latex?i
            "i")-th window and
            ![f\_2](https://latex.codecogs.com/png.latex?f_2 "f_2") be
            the pdf of (i + 1)-th window
        2.  K-L divergence
            ![d\_i](https://latex.codecogs.com/png.latex?d_i "d_i") from
            ![f\_2](https://latex.codecogs.com/png.latex?f_2 "f_2") to
            ![f\_1](https://latex.codecogs.com/png.latex?f_1 "f_1")
        3.  Set ![D =
            d\_i](https://latex.codecogs.com/png.latex?D%20%3D%20d_i
            "D = d_i").
        4.  ![f^{\\prime} =
            f\_1](https://latex.codecogs.com/png.latex?f%5E%7B%5Cprime%7D%20%3D%20f_1
            "f^{\\prime} = f_1")
    2.  else
        1.  ![f^{\\prime}](https://latex.codecogs.com/png.latex?f%5E%7B%5Cprime%7D
            "f^{\\prime}") be the pdf of normal and
            ![f\_2](https://latex.codecogs.com/png.latex?f_2 "f_2") be
            the pdf of (i + 1)-th window
        2.  K-L divergence
            ![d\_i](https://latex.codecogs.com/png.latex?d_i "d_i") from
            ![f\_2](https://latex.codecogs.com/png.latex?f_2 "f_2") to
            ![f^{\\prime}](https://latex.codecogs.com/png.latex?f%5E%7B%5Cprime%7D
            "f^{\\prime}")
        3.  Set ![D =
            d\_i](https://latex.codecogs.com/png.latex?D%20%3D%20d_i
            "D = d_i").
2.  output: ![\\{ d\_1, \\ldots, d\_{w - 1}
    \\}](https://latex.codecogs.com/png.latex?%5C%7B%20d_1%2C%20%5Cldots%2C%20d_%7Bw%20-%201%7D%20%5C%7D
    "\\{ d_1, \\ldots, d_{w - 1} \\}")

### Dynamic lambda algorithm

  - Data: univariate series of size
    ![n](https://latex.codecogs.com/png.latex?n "n")
  - Input: window size ![win](https://latex.codecogs.com/png.latex?win
    "win"), jump size for sliding window
    ![jump](https://latex.codecogs.com/png.latex?jump "jump"), threshold
    increment
    ![\\lambda^{\\prime}](https://latex.codecogs.com/png.latex?%5Clambda%5E%7B%5Cprime%7D
    "\\lambda^{\\prime}"), threshold updating constant
    ![\\epsilon](https://latex.codecogs.com/png.latex?%5Cepsilon
    "\\epsilon")

The threshold is initialized by

  
![\\lambda = \\lambda^{\\prime}
\\epsilon](https://latex.codecogs.com/png.latex?%5Clambda%20%3D%20%5Clambda%5E%7B%5Cprime%7D%20%5Cepsilon
"\\lambda = \\lambda^{\\prime} \\epsilon")  

If ![D \<
\\lambda](https://latex.codecogs.com/png.latex?D%20%3C%20%5Clambda
"D \< \\lambda"), it is updated by

  
![\\lambda = \\lambda^{\\prime} (d\_{j - 2} +
\\epsilon)](https://latex.codecogs.com/png.latex?%5Clambda%20%3D%20%5Clambda%5E%7B%5Cprime%7D%20%28d_%7Bj%20-%202%7D%20%2B%20%5Cepsilon%29
"\\lambda = \\lambda^{\\prime} (d_{j - 2} + \\epsilon)")  

1.  ![\\lambda = \\lambda^{\\prime}
    \\epsilon](https://latex.codecogs.com/png.latex?%5Clambda%20%3D%20%5Clambda%5E%7B%5Cprime%7D%20%5Cepsilon
    "\\lambda = \\lambda^{\\prime} \\epsilon")
2.  ![d\_1](https://latex.codecogs.com/png.latex?d_1 "d_1") K-L from
    ![1](https://latex.codecogs.com/png.latex?1 "1")-st to
    ![2](https://latex.codecogs.com/png.latex?2 "2")-nd window
3.  ![d\_2](https://latex.codecogs.com/png.latex?d_2 "d_2") K-L from
    ![2](https://latex.codecogs.com/png.latex?2 "2")-nd to
    ![3](https://latex.codecogs.com/png.latex?3 "3")-rd window
4.  For ![i
    \\leftarrow 3](https://latex.codecogs.com/png.latex?i%20%5Cleftarrow%203
    "i \\leftarrow 3") to (w - 1) do
    1.  if ![D \<
        \\lambda](https://latex.codecogs.com/png.latex?D%20%3C%20%5Clambda
        "D \< \\lambda") then
        1.  ![f\_1](https://latex.codecogs.com/png.latex?f_1 "f_1") be
            the pdf of ![i](https://latex.codecogs.com/png.latex?i
            "i")-th window and
            ![f\_2](https://latex.codecogs.com/png.latex?f_2 "f_2") be
            the pdf of (i + 1)-th window
        2.  K-L divergence
            ![d\_i](https://latex.codecogs.com/png.latex?d_i "d_i") from
            ![f\_2](https://latex.codecogs.com/png.latex?f_2 "f_2") to
            ![f\_1](https://latex.codecogs.com/png.latex?f_1 "f_1")
        3.  ![\\lambda = \\lambda^{\\prime} (d\_{j - 2} +
            \\epsilon)](https://latex.codecogs.com/png.latex?%5Clambda%20%3D%20%5Clambda%5E%7B%5Cprime%7D%20%28d_%7Bj%20-%202%7D%20%2B%20%5Cepsilon%29
            "\\lambda = \\lambda^{\\prime} (d_{j - 2} + \\epsilon)")
        4.  Set ![D =
            d\_i](https://latex.codecogs.com/png.latex?D%20%3D%20d_i
            "D = d_i").
        5.  ![f^{\\prime} =
            f\_1](https://latex.codecogs.com/png.latex?f%5E%7B%5Cprime%7D%20%3D%20f_1
            "f^{\\prime} = f_1")
    2.  else
        1.  ![f^{\\prime}](https://latex.codecogs.com/png.latex?f%5E%7B%5Cprime%7D
            "f^{\\prime}") be the pdf of normal and
            ![f\_2](https://latex.codecogs.com/png.latex?f_2 "f_2") be
            the pdf of (i + 1)-th window
        2.  K-L divergence
            ![d\_i](https://latex.codecogs.com/png.latex?d_i "d_i") from
            ![f\_2](https://latex.codecogs.com/png.latex?f_2 "f_2") to
            ![f^{\\prime}](https://latex.codecogs.com/png.latex?f%5E%7B%5Cprime%7D
            "f^{\\prime}")
        3.  Set ![D =
            d\_i](https://latex.codecogs.com/png.latex?D%20%3D%20d_i
            "D = d_i").
5.  output: ![\\{ d\_1, \\ldots, d\_{w - 1}
    \\}](https://latex.codecogs.com/png.latex?%5C%7B%20d_1%2C%20%5Cldots%2C%20d_%7Bw%20-%201%7D%20%5C%7D
    "\\{ d_1, \\ldots, d_{w - 1} \\}")

-----

# References

<div id="refs" class="references">

<div id="ref-Cho:2019ji">

Cho, Jinwoo, Shahroz Tariq, Sangyup Lee, Young Geun Kim, and Simon Woo.
2019. “Contextual Anomaly Detection by Correlated Probability
Distributions using Kullback-Leibler Divergence.” *WORKSHOP ON MINING
AND LEARNING FROM TIME SERIES*.

</div>

<div id="ref-Yun:2018di">

Yun, Jeong-Han, Yoonho Hwang, Woomyo Lee, Hee-Kap Ahn, and Sin-Kyu Kim.
2018. “Statistical Similarity of Critical Infrastructure Network Traffic
Based on Nearest Neighbor Distances.” In *Research in Attacks,
Intrusions, and Defenses*, 1–23. Cham: Springer International
Publishing.

</div>

</div>
