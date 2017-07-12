# Truncated Distributions
`(truncated d a b)` returns a distribution d truncated to the interval (a,b], which behaves like any other distribution.


```
(truncated (normal 0 4) 1 2)
#distributions.core.TruncatedDistribution{:distribution #distributions.core.Normal{:mean 0, :variance 4}, :lower 1, :upper 2, :F-lower 0.691462461274013, :F-upper 0.841344746068543}
```
Notice that the truncated distribution is implemented as a record containing the original distribution, the lower and upper truncation limits, and the original cdf evaluated thereat as these values are frequently used in the pdf cdf and icdf functions.

```
λ> (cdf (truncated (normal 0 4) 1 2) 1)
0.0

λ> (cdf (truncated (normal 0 3) 1 2) 1.4393822324537189)
0.49999999999999967

λ> (cdf (truncated (normal 0 4) 1 2) 2)
1.0

λ> (icdf (truncated (normal 0 3) 1 2) 0)
0.9999999999999997

λ> (icdf (truncated (normal 0 3) 1 2) 0.5)
1.4393822324537189

λ> (icdf (truncated (normal 0 3) 1 2) 1)
2.0000000000000004

```

Sampling is achcieved through the inverse-cdf method.
```
λ> (sample (truncated (gamma 4 3) 5 6))
5.424741741923803

λ> (sample (truncated (gamma 4 3) 5 6) 4)
(5.0826669400151445 5.957356728304568 5.393586364632808 5.2298513236000135)
```

pdf and log-pdf are as follows
```
λ> (pdf (truncated (normal 0 4) 1 2) 1)
0

λ> (pdf (truncated (normal 0 4) 1 2) 1.5)
1.0045798026352024

λ> (pdf (truncated (normal 0 4) 1 2) 2)
0.8072025484894866

λ> (log-pdf (truncated (normal 0 4) 1 2) 1)
-Infinity

λ> (log-pdf (truncated (normal 0 4) 1 2) 2)
-0.21418065275063713
```
!!! note "pdf evaluated at lower limit"
    Note that because the truncated interval is (a,b], the pdf evaluated at a is zero.


The support is now determined by the truncation interval
```
λ> (support-lower (truncated (gamma 1 2) 4 5))
4

λ> (support-upper (truncated (gamma 1 2) 4 5))
5
```
