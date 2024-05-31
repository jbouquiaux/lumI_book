---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Theory

```{code-cell}
from lineshape_slider import get_plotly_lineshape
%matplotlib inline
```
```{code-cell}
get_plotly_lineshape(slider="frequency",num_step=10)
```

```{code-cell}
get_plotly_lineshape(slider="Delta_Q",num_step=10)
```



