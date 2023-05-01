# Polylogs in Python!
Polylogarithms in Numpy/Scipy.

A collection of functions used to calculate polylogarithms in Python using Numpy & Scipy.

```python
# Real only Polylog
result = polylog_real(n, z)

# Complex Polylog
result = polylog(n, z)

# Recursive Polylog
result = polylog_rec(n, z)


# All functions also have a precision input, dictating the number of decimal points evaluation should be correct to.
```

This allows the evaluation of complex functions such as:

$$\int_0^k \frac{u^3}{e^u-1} \mathrm{d}x = \frac{\pi^4}{15} + k^3\log({1-e^k}) - 3k^2 \mathrm{Li}_2(e^{-k})-6k\mathrm{Li}_3(e^{-k})-6\mathrm{Li}_4(e^{-k})$$

where $\mathrm{Li}_n(z)$ is the polylogarithm.

These functions often appear in radiometrical calculations and across physics & engineering.