# Notebook Review

This note summarizes the strongest issues to fix or explain before turning the notebook into a report.

## Key Findings

1. The stated linking-length search range does not match the implemented range.
The markdown says the search is over `5 Mpc <= r_c <= 25 Mpc`, but the code actually uses `np.arange(1, 11, .5)`.
Reference: [ACN - from scratch copy.txt](/Users/ali/Course%20Materials/Complex%20networks/PROJECT/ACN%20-%20from%20scratch%20copy.txt:188) and [ACN - from scratch copy.txt](/Users/ali/Course%20Materials/Complex%20networks/PROJECT/ACN%20-%20from%20scratch%20copy.txt:218)
Impact: the reported `r_c` is only valid inside the narrower scan window, so the report must either fix the code or describe the real search interval.

2. The percolation-threshold rule is a heuristic, not a validated threshold estimate.
The code chooses `r_c` from the maximum numerical gradient of GCC coverage. That is a reasonable first pass, but it can be sensitive to step size, noise, and the chosen scan interval.
Reference: [ACN - from scratch copy.txt](/Users/ali/Course%20Materials/Complex%20networks/PROJECT/ACN%20-%20from%20scratch%20copy.txt:241)
Impact: in the report, describe this as a heuristic estimate unless you add a stronger selection rule and sensitivity test.

3. The clustering validation is asymmetric.
The code removes field galaxies only from the true labels (`True_Cluster_ID != 0`) but keeps all predicted Louvain labels, including singleton or tiny communities.
Reference: [ACN - from scratch copy.txt](/Users/ali/Course%20Materials/Complex%20networks/PROJECT/ACN%20-%20from%20scratch%20copy.txt:357)
Impact: the NMI result is useful, but it should be described carefully. A cleaner comparison would also filter predicted groups by a minimum size and explicitly document which galaxies are included in the score.

4. The modularity value is probably not enough by itself to claim excellent physical clustering.
For a sparse graph with many disconnected or weakly connected components, modularity can become very high even when the partition is not physically meaningful.
Reference: [ACN - from scratch copy.txt](/Users/ali/Course%20Materials/Complex%20networks/PROJECT/ACN%20-%20from%20scratch%20copy.txt:333)
Impact: keep the modularity number, but pair it with NMI, component statistics, and a limitation note.

5. The projection plot asks for a legend without labeled artists.
Matplotlib already warns about this in the saved output.
Reference: [ACN - from scratch copy.txt](/Users/ali/Course%20Materials/Complex%20networks/PROJECT/ACN%20-%20from%20scratch%20copy.txt:1908)
Impact: either remove the legend call or assign labels to the plotted clusters.

6. The degree-distribution interpretation is too strong as written.
The notebook says a straight descending line on a log-log plot means a true cosmic web or scale-free structure. That conclusion needs a fitted model and goodness-of-fit test; a visual trend alone is not enough.
Reference: [ACN - from scratch copy.txt](/Users/ali/Course%20Materials/Complex%20networks/PROJECT/ACN%20-%20from%20scratch%20copy.txt:2062)
Impact: phrase this section as exploratory topological evidence, not proof of a power law.

## What Is Already Working Well

- The data-preparation flow is coherent: fetch catalog, apply redshift and magnitude cuts, and carry labels through the pipeline.
- The coordinate conversion uses `astropy`, which is a solid scientific choice.
- The KDTree-based edge search is appropriate for a large galaxy sample.
- The notebook already produces useful report figures and concrete summary numbers.

## Recommended Next Steps

1. Decide on the real `r_c` search interval and make the markdown and code consistent.
2. Reframe the `r_c` choice as a heuristic unless you add a sensitivity analysis.
3. Add one small validation table:
`nodes`, `edges`, `number of communities`, `largest component fraction`, `modularity`, `NMI`.
4. Tone down the degree-distribution claim unless you run a proper power-law fit.
5. Clean the figure labels and captions so they are ready to paste into the report.

## Suggested Report Positioning

The current notebook is already strong enough for a course report if you present it honestly as:

- a network-based reconstruction of galaxy structure in a magnitude-limited SDSS slice,
- using a KDTree friends-of-friends style graph,
- with Louvain communities as a topology-based proxy for superclusters,
- and with NMI used as an external validation against catalog labels.

That framing is both defensible and aligned with what the code actually does.
