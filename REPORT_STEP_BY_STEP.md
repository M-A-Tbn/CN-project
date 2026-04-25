# Report Guide

Use this as the writing scaffold for the project report. The text below is intentionally half-outline and half-draft so you can expand it section by section.

## Working Title

Network-Based Reconstruction of Large-Scale Galaxy Structure in an SDSS DR12 Slice

## 1. Abstract

What to include:

- the scientific goal,
- the dataset,
- the main method,
- the main result,
- one sentence on validation.

Starter draft:

This project studies large-scale structure in a magnitude-limited slice of the SDSS DR12 galaxy catalog using a network-based approach. After selecting galaxies in the redshift range `0.04 <= z <= 0.10`, I converted sky coordinates into 3D comoving Cartesian coordinates and constructed a graph in which nearby galaxies were connected within a linking length `r_c`. I estimated `r_c` from the growth of the giant connected component, then applied Louvain community detection to identify topology-based superclusters. The resulting network contained about `145,583` galaxies and `602,240` edges, and the detected communities achieved a normalized mutual information of about `0.805` against catalog cluster labels. These results suggest that graph methods can recover meaningful large-scale structure, while also highlighting sensitivity to the choice of linking length and validation procedure.

## 2. Introduction

What to include:

- why the cosmic web matters,
- why networks are useful,
- what question you are trying to answer.

Starter draft:

The large-scale distribution of galaxies is not random. Galaxies form filaments, clusters, walls, and voids that together make up the cosmic web. Because this structure is fundamentally relational, network theory provides a natural framework for analysis: galaxies can be treated as nodes, and physically nearby galaxies can be linked by edges. The main goal of this project is to test whether a simple proximity-based network, built from SDSS DR12 data, can recover meaningful galaxy group structure and produce interpretable large-scale patterns.

## 3. Data

What to include:

- source catalog,
- why that catalog is suitable,
- the selected columns,
- the cuts you applied.

Starter draft:

The data were obtained from the Tempel SDSS DR12 galaxy-group catalog through Vizier. From the original table, I retained the object identifier, right ascension, declination, observed redshift, apparent `r`-band magnitude, and catalog group identifier. To work with a controlled volume and avoid severe selection bias, I restricted the sample to galaxies with redshifts between `0.04` and `0.10`. I also imposed a magnitude-limited cut by computing the absolute-magnitude threshold corresponding to the survey limit at the far edge of the volume.

## 4. Preprocessing and Coordinate Transformation

What to include:

- cosmology used,
- distance modulus and magnitude limit idea,
- conversion from spherical to Cartesian coordinates.

Starter draft:

I adopted a flat Lambda-CDM cosmology with `H0 = 70 km s^-1 Mpc^-1` and `Om0 = 0.3`. Using the SDSS apparent-magnitude limit `m = 17.77`, I computed the distance modulus at `z = 0.10` and derived an absolute-magnitude threshold of about `M = -20.55`. Only galaxies brighter than this limit were kept, ensuring that the final sample was approximately volume-limited across the selected redshift interval. I then converted the observed sky positions and comoving distances into 3D Cartesian coordinates `(X, Y, Z)` using `astropy.coordinates.SkyCoord`.

## 5. Network Construction

What to include:

- node definition,
- edge definition,
- KDTree use,
- how `r_c` was chosen.

Starter draft:

Each galaxy was represented as a node in a spatial graph. Edges were added between galaxy pairs separated by less than a linking length `r_c`, which was found efficiently with a KDTree. To estimate a suitable `r_c`, I measured the fraction of galaxies belonging to the giant connected component as a function of trial linking length and selected the value where this curve rose most sharply. In the current notebook, this produces an `r_c` close to `5.1 Mpc`, although this value should be interpreted as a heuristic estimate because it depends on the scan range and step size.

Important note for the final version:

- Make the search interval in the text match the code.
- If you keep the current method, explicitly call it a heuristic percolation-threshold estimate.

## 6. Community Detection

What to include:

- why Louvain,
- what modularity means,
- what the communities represent physically.

Starter draft:

After fixing the linking length, I built the final galaxy graph and applied the Louvain algorithm to detect communities. Louvain partitions the graph by maximizing modularity, so the detected groups are regions with dense internal connections and weaker external links. In a cosmological context, these communities can be interpreted as topology-based superclusters or dense substructures within the cosmic web. For the current network, Louvain identified about `16,375` communities with a modularity of about `0.994`.

## 7. Validation

What to include:

- comparison to Tempel labels,
- what NMI measures,
- what the current score means,
- one sentence on limitations.

Starter draft:

To evaluate the detected communities, I compared the Louvain labels to the group identifiers from the Tempel catalog using normalized mutual information (NMI). NMI measures how much information is shared between two partitions, with higher values indicating stronger agreement. The current notebook gives an NMI of about `0.805`, which suggests substantial overlap between the network-based partition and the catalog labels. However, this comparison should be interpreted carefully because the scoring currently excludes field galaxies only on the catalog side and therefore is not a perfectly symmetric validation.

## 8. Visual Results

What to include:

- 3D scatter,
- 2D projections,
- filament plot,
- what each figure is meant to show.

Starter draft:

I used three types of visualization to interpret the network. First, an interactive 3D scatter plot highlights the largest detected communities in full spatial context. Second, 2D `X-Y` and `X-Z` projections reveal the anisotropic geometry of the detected structures while preserving the low-density background field. Third, a filament-style plot overlays sampled network edges on the galaxy distribution, making the web-like connectivity visible in projection. Together, these figures show that the graph captures both dense nodes and extended connecting structures.

## 9. Topological Analysis

What to include:

- degree,
- hubs,
- degree-distribution plot,
- cautious interpretation.

Starter draft:

To characterize the network topology, I computed the degree of every galaxy, where degree is the number of neighboring galaxies connected to that node. The highest-degree galaxies act as hubs and are likely associated with dense cluster cores. I also plotted the degree distribution on log-log axes to inspect whether the network shows heterogeneous connectivity. This plot is useful as an exploratory diagnostic, but it should not be presented as proof of scale-free behavior without a dedicated statistical fit.

## 10. Limitations

Good points to include:

- sensitivity to linking length,
- redshift-space distortions,
- magnitude-limited selection,
- community detection is algorithm-dependent,
- validation is not perfectly symmetric.

Starter draft:

This analysis has several limitations. The graph depends strongly on the chosen linking length, and the current threshold is selected by a heuristic rule rather than a formal optimization. The sample is only approximately volume-limited and is still affected by observational selection. In addition, observed redshifts include peculiar-velocity distortions, so the 3D geometry is not purely real-space structure. Finally, Louvain communities are algorithmic partitions and should be interpreted as network-based proxies for physical superclusters rather than exact ground truth.

## 11. Conclusion

Starter draft:

This project shows that a relatively simple graph-construction pipeline can recover meaningful large-scale galaxy structure from SDSS data. By combining cosmological preprocessing, KDTree-based neighbor search, and Louvain community detection, the notebook identifies coherent communities that agree reasonably well with catalog labels. The strongest next improvement would be a more careful treatment of the linking-length selection and a cleaner validation protocol. Even in its current form, the analysis demonstrates that network science provides a useful language for studying the cosmic web.

## 12. Figures to Include

Recommended figure order:

1. GCC fraction versus linking length.
2. 3D map of the largest detected communities.
3. `X-Y` and `X-Z` projections.
4. Filament-style projected network figure.
5. Degree distribution.

## 13. Small Table to Add

Add a one-table summary with:

- number of galaxies,
- redshift range,
- magnitude limit,
- chosen `r_c`,
- number of edges,
- number of Louvain communities,
- modularity,
- NMI.

## 14. Best Next Writing Order

Write the report in this order:

1. Data
2. Methods
3. Validation
4. Results
5. Introduction
6. Conclusion
7. Abstract

That order is easier because the technical content is already in the notebook.
