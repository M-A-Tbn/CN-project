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

## How To Make This Stronger Than Riazi's Report

To make your project stronger than hers across the board, you do not need to copy her style. You need to make your pipeline more defensible, better validated, and more tightly argued.

Right now, your project is already more ambitious technically. The places where her work is stronger are mostly about discipline: the question is sharper, the method matches the question closely, and the validation is easier to trust. So the path forward is to keep your richer pipeline and tighten the weak links.

### What to improve first

`1. Make your research question brutally clear`

Her report wins because it has one crisp claim: communities correspond to subfields.
You should define one equally sharp claim, for example:

- a proximity-based galaxy network can recover known cluster structure from SDSS data
- community detection on a spatial galaxy graph identifies meaningful large-scale structure
- graph topology captures physically relevant cosmic organization beyond simple spatial clustering

Pick one and make every section serve that claim.

`2. Strengthen the linking-length choice`

This is your biggest methodological weakness right now. Your whole graph depends on `r_c`, and the current choice is heuristic.

To beat her here, do this:

- run a sensitivity analysis over several `r_c` values
- show how `NMI`, number of communities, giant component size, and modularity change with `r_c`
- justify your final `r_c` using both percolation behavior and validation performance

If you do that, your graph construction becomes more rigorous than her ready-made citation graph workflow.

`3. Clean up the validation`

Your validation idea is actually stronger than hers, but it needs to be executed more carefully.

Improve it by:

- making the catalog-vs-predicted comparison symmetric
- clearly defining which galaxies are included in the score
- reporting more than one metric: `NMI`, `ARI`, purity, maybe size-matched overlap
- comparing against a baseline method, not just Louvain alone

If you do that, your validation becomes clearly better than her Fisher test.

`4. Compare methods, not just run one`

She uses one method and explains why it fits.
You can go further by comparing alternatives.

For example:

- Louvain vs Leiden vs connected components vs DBSCAN / FoF-style grouping
- different linking lengths
- maybe weighted vs unweighted graphs

Then report:

- which method best matches catalog structure
- which method gives the most stable communities
- which method is most physically interpretable

That would make your methods section much stronger than hers.

`5. Be more careful with interpretation`

A place where her report is stronger is restraint. She usually says only what the evidence supports.

To improve yours:

- do not call the degree plot evidence of a power law unless you fit it properly
- do not treat high modularity as proof of physical truth
- clearly separate exploratory observations from validated conclusions

Paradoxically, sounding a little more cautious will make the project feel more scientific and therefore stronger.

### How to make each field better

`Data`

You are already strong here. To go further:

- explain why your redshift slice is chosen
- justify the magnitude cut more explicitly
- mention selection effects and redshift-space distortions

This can become better than hers without much extra work.

`Methods`

This is where you can win decisively if you:

- justify `r_c`
- compare community methods
- explain why Louvain is appropriate for your graph
- discuss limitations of graph construction

`Validation`

This is the highest upside section. Add:

- symmetric filtering
- multiple validation metrics
- baseline comparison
- sensitivity across parameters

`Results`

Keep your rich figures, but organize them around the main claim:

- does the method recover real structure?
- where does it succeed?
- where does it fail?

That will make your results more convincing than just visually impressive.

`Discussion`

You can outperform her by explicitly discussing:

- what graph-based recovery captures well
- where it diverges from catalog truth
- whether communities correspond to clusters, superclusters, or something in between

That kind of conceptual clarity will lift the whole project.

### The shortest honest answer

If you want your work to be better than hers in all areas, focus on these four upgrades:

1. Sharpen the main question.
2. Justify `r_c` with a sensitivity study.
3. Make validation cleaner and more extensive.
4. Interpret results more cautiously and precisely.

You already have more technical depth. What you need now is stronger scientific control.
