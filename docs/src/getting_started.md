# Getting Started

This guide introduces the basic workflow of **VisGraphs.jl** for constructing and analyzing visibility graphs from time series data.

We cover:

* Horizontal Visibility Graph (HVG)
* Natural Visibility Graph (NVG)
* Weighted variants
* Basic graph analysis tools

---

## Installation

If you haven’t installed **VisGraphs.jl** yet:

```julia
using Pkg
Pkg.add(url="https://github.com/Mangojoghurt/VisGraphs.jl")
```

Then load the package:

```@example main
using VisGraphs
```

---

## Example Time Series

We start by generating a simple noisy oscillatory signal using the built-in utility function from **VisGraphs.jl**.

```@example main
x = generate_noisy_sine(50, 0.1)
```

This produces a sine wave with additive Gaussian noise, which is useful for demonstrating visibility graph constructions on realistic, non-smooth signals.

---

## Horizontal Visibility Graph (HVG)

The **Horizontal Visibility Graph (HVG)** connects points that can “see” each other using a horizontal line-of-sight criterion.

### Construct the HVG

```@example main
edges_hvg = hvg(x)
```

### Plot the HVG

```@example main
plot_hvg(x)
```

---

## Natural Visibility Graph (NVG)

The **Natural Visibility Graph (NVG)** uses a geometric criterion based on straight-line visibility between points.

### Construct the NVG

```@example main
edges_nvg = nvg(x)
```

### Plot the NVG

```@example main
plot_nvg(x)
```

---

## Weighted Visibility Graphs

Weighted variants attach one number to each edge: the angle, in radians,
of the line connecting the two nodes, computed as
`atan(x[j] - x[i], j - i)`. This combines the vertical difference in
value and the horizontal distance in time/index into a single quantity
bounded within `(-π/2, π/2)` — a steep, fast change between two points
produces an angle with large magnitude, while a gradual, slow change
produces an angle close to zero.

```@example main
edges_whvg = whvg(x)
```
```@example main
edges_wnvg = wnvg(x)
```

```@example main
plot_whvg(x)
```
```@example main
plot_wnvg(x)
```

---

## Graph Analysis

Once constructed, visibility graphs can be analyzed using standard graph-theoretic tools.

### Adjacency Matrix

```@example main
A = adjacency_matrix(edges_hvg, length(x))
```

### Degree Distribution

```@example main
degrees, dist = degree_distribution(edges_hvg, length(x))
```

### Laplacian Matrix

```@example main
L = laplacian_matrix(edges_hvg, length(x))
```

### Clustering Coefficient

```@example main
local_c, global_c = clustering_coefficient(edges_hvg, length(x))
```

### Average Path Length

```@example main
average_path_length(edges_hvg, length(x))
```

!!! note
    `adjacency_matrix` and `laplacian_matrix` return **sparse** matrices
    (`SparseMatrixCSC`), since visibility graphs typically have far fewer
    edges than a fully dense graph would. All standard matrix operations
    (indexing, `size`, `sum`, `==`, etc.) work exactly as they would on a
    dense matrix.

---

## Comparing HVG and NVG

It is often useful to compare the structural differences between HVG and NVG representations.

```@example main
plot_hvg(x)
```
```@example main
plot_nvg(x)
```

---

## Next Steps

You can now:

* Try different time series (chaotic, financial, periodic)
* Compare HVG vs NVG structures
* Explore spectral properties via the Laplacian
* Use weighted graphs for richer feature extraction
