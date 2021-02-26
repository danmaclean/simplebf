`simplebf` is a package of convenience functions for Bayes Factor
analysis of two samples. The package is developed as a teaching aid for
those new to Bayesian Inference and working with R.

The functions are thin wrappers around functions in `BayesFactor` on
CRAN [here](https://cran.r-project.org/package=BayesFactor) and project
webpage <http://bayesfactor.blogspot.com/>. Please use that in
preference.

Using `simplebf` for a *t*-test
-------------------------------

To carry out Bayes Factor analysis on two named groups in a tidy
dataframe use `named_pair_ttestbf()`, pass it the column names for the
group and data measurements, the name of the group to be considered as
the test group, and the name of the group to be considere the control
group.

``` r
library(simplebf)
result <- named_pair_ttestbf(PlantGrowth, 
                             group_col = "group", data_col = "weight",
                             control = "ctrl", test = "trt2"
                             )
```

The function returns a one row dataframe, the columns are as follows:

    1. `control_group` - the name of the group considered control
    2. `test_group` - the group considered test
    3. `h_0` - statement of the $H_0$ hypothesis
    4. `h_1` - statement of the $H_1$ hypothesis
    5. `BayesFactor` - the resulting Bayes Factor from $\frac{H_{1}/H_{0}}$
    6. `odds_h_1` - statement of the Bayes Factor in terms of odds $H_0:H_1$
    7. `summary` - summary statement of relative evidence for the hypothesis

``` r
result %>% knitr::kable()
```

<table>
<colgroup>
<col style="width: 9%" />
<col style="width: 7%" />
<col style="width: 13%" />
<col style="width: 16%" />
<col style="width: 8%" />
<col style="width: 13%" />
<col style="width: 31%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">control_group</th>
<th style="text-align: left;">test_group</th>
<th style="text-align: left;">h_0</th>
<th style="text-align: left;">h_1</th>
<th style="text-align: right;">BayesFactor</th>
<th style="text-align: left;">odds_h_1</th>
<th style="text-align: left;">summary</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">ctrl</td>
<td style="text-align: left;">trt2</td>
<td style="text-align: left;">trt2 equal to ctrl</td>
<td style="text-align: left;">trt2 greater than ctrl</td>
<td style="text-align: right;">3.387166</td>
<td style="text-align: left;">1:3.38716599374275</td>
<td style="text-align: left;">Substantial evidence for H_1 compared to H_0</td>
</tr>
</tbody>
</table>

### Changing the scale of the prior and the nature *H*<sub>1</sub>

Only one option for the `BayesFactor::ttestBF()` function can be
changed. The value of `rscale` for the prior distribution can be set to
one of `medium`, `wide` or `ultrawide`. The hypothesis to test as
*H*<sub>1</sub> can be set to one of `test_greater_than_control`,
`test_smaller_than_control`,`test_not_equal_to_control`.

``` r
result <- named_pair_ttestbf(PlantGrowth, 
                             group_col = "group", data_col = "weight",
                             control = "ctrl", test = "trt2",
                             rscale = "medium",
                             h_1 = "test_greater_than_control"
                             )
```

Comparing all groups in a dataframe
-----------------------------------

All levels of the group column can be compared pairwise using the
`allpairs_ttestbf()` function. This is called similarly to the
`named_pair_ttestbf()` but without the need to specify control and test
groups as all pairs are done automatically. A dataframe with the same
columns as above and one row for each pair (in `A:B` *and* `B:A` order)
is returned. `rscale` and `h_1` can be used here too.

``` r
result <- allpairs_ttestbf(PlantGrowth, group_col = "group", data_col = "weight")
result %>% knitr::kable()
```

<table>
<colgroup>
<col style="width: 9%" />
<col style="width: 7%" />
<col style="width: 13%" />
<col style="width: 15%" />
<col style="width: 8%" />
<col style="width: 13%" />
<col style="width: 31%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">control_group</th>
<th style="text-align: left;">test_group</th>
<th style="text-align: left;">h_0</th>
<th style="text-align: left;">h_1</th>
<th style="text-align: right;">BayesFactor</th>
<th style="text-align: left;">odds_h_1</th>
<th style="text-align: left;">summary</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">trt1</td>
<td style="text-align: left;">ctrl</td>
<td style="text-align: left;">ctrl equal to trt1</td>
<td style="text-align: left;">ctrl greater than trt1</td>
<td style="text-align: right;">1.0833663</td>
<td style="text-align: left;">1:1.08336631965942</td>
<td style="text-align: left;">Anecdotal evidence for H_1 compared to H_0</td>
</tr>
<tr class="even">
<td style="text-align: left;">trt2</td>
<td style="text-align: left;">ctrl</td>
<td style="text-align: left;">ctrl equal to trt2</td>
<td style="text-align: left;">ctrl greater than trt2</td>
<td style="text-align: right;">0.1622109</td>
<td style="text-align: left;">1:0.162210898064748</td>
<td style="text-align: left;">Substantial evidence for H_0 compared to H_1</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ctrl</td>
<td style="text-align: left;">trt1</td>
<td style="text-align: left;">trt1 equal to ctrl</td>
<td style="text-align: left;">trt1 greater than ctrl</td>
<td style="text-align: right;">0.2167156</td>
<td style="text-align: left;">1:0.216715621920709</td>
<td style="text-align: left;">Substantial evidence for H_0 compared to H_1</td>
</tr>
<tr class="even">
<td style="text-align: left;">trt2</td>
<td style="text-align: left;">trt1</td>
<td style="text-align: left;">trt1 equal to trt2</td>
<td style="text-align: left;">trt1 greater than trt2</td>
<td style="text-align: right;">0.1363151</td>
<td style="text-align: left;">1:0.13631506044811</td>
<td style="text-align: left;">Substantial evidence for H_0 compared to H_1</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ctrl</td>
<td style="text-align: left;">trt2</td>
<td style="text-align: left;">trt2 equal to ctrl</td>
<td style="text-align: left;">trt2 greater than ctrl</td>
<td style="text-align: right;">3.3871660</td>
<td style="text-align: left;">1:3.38716599374275</td>
<td style="text-align: left;">Substantial evidence for H_1 compared to H_0</td>
</tr>
<tr class="even">
<td style="text-align: left;">trt1</td>
<td style="text-align: left;">trt2</td>
<td style="text-align: left;">trt2 equal to trt1</td>
<td style="text-align: left;">trt2 greater than trt1</td>
<td style="text-align: right;">12.6444900</td>
<td style="text-align: left;">1:12.6444900067635</td>
<td style="text-align: left;">Strong evidence for H_1 compared to H_0</td>
</tr>
</tbody>
</table>
