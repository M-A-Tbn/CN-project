# Network-Based Reconstruction of Galaxy Structure

## Phase 1: Importing data from SDSS DR12 database.

```python
# importing the labeled Galaxy Clusters via Vizier from Tempel's SDSS DR12 catalog
from astroquery.vizier import Vizier
import pandas as pd
vizier = Vizier(row_limit=-1)
print("Querying Vizier for galaxy cluster data...")
catalogs = vizier.get_catalogs("J/A+A/602/A100")
print("Found these tables: \n")
print(catalogs.keys())

#defines the cosmology
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3) 

# SDSS apparent magnitude limit for the main galaxy sample
apparent_magnitude_limit = 17.77

# Since we can't work with the entire database, we choose a limited volume of the space to work with.(Redshifts between 0.04 and 0.1)
z_min = 0.04
z_max = 0.1

# then we must Calculate the absolute magnitude limit at the FAR edge of our limited volume.
# Any galaxy fainter than this wouldn't be detectable at z_max. (M = m - mu)
mu_max = cosmo.distmod(z_max).value
M_limit = apparent_magnitude_limit - mu_max
print(f"Absolute magnitude limit at z={z_max}: {M_limit:.2f}")

# 1. Grab the first table (the galaxy list)
df = catalogs[0].to_pandas()

# 2. Keep only the columns we actually need for the math and the network
cosmic_df = df[['objID', 'RAJ2000', 'DEJ2000', 'zobs', 'rmag', 'GroupID']].copy()

# 3. Rename them to simple, standard variables
cosmic_df.rename(columns={
    'RAJ2000': 'ra', 
    'DEJ2000': 'dec', 
    'zobs': 'z', 
    'GroupID': 'True_Cluster_ID',
    'rmag': 'm'
}, inplace=True)

# calculate the absolute magnitude for each galaxy (M = m - mu)
cosmic_df['M'] = cosmic_df['m'] - cosmo.distmod(cosmic_df['z']).value

# 4. Filter the data to a specific "slice" of the universe (Redshift between 0.04 and 0.1)
# also filter to only include galaxies brighter than the derived absolute magnitude limit

my_galaxies = cosmic_df[
    (cosmic_df['z'] >= z_min) & 
    (cosmic_df['z'] <= z_max) &
    (cosmic_df['M'] <= M_limit)
    ].copy()

# Reset the index to be clean
my_galaxies.reset_index(drop=True, inplace=True)

print("Data prepped! Here are the first 5 galaxies:")
print(my_galaxies.head())
```

## Phase 2: Coordinate Transformation
At this point, We need to convert our spherical coordinates (d, $\alpha, \delta$) to cartesian coordinates (X,Y,Z).

```python
from astropy import units as u
from astropy.coordinates import SkyCoord

#calculate distance from redshift
my_galaxies['distance_mpc'] = cosmo.comoving_distance(my_galaxies['z']).to(u.Mpc).value

#summon the skycoord and let it do its magic
coords = SkyCoord(
    ra=my_galaxies['ra'].values * u.degree,
    dec=my_galaxies['dec'].values * u.degree,
    distance=my_galaxies['distance_mpc'].values * u.Mpc
)

#extract the 3D coordinates
my_galaxies['X'] = coords.cartesian.x.value
my_galaxies['Y'] = coords.cartesian.y.value
my_galaxies['Z'] = coords.cartesian.z.value

print("Phase 2 Complete! Here are the 3D coordinates:")
print(my_galaxies[['objID', 'X', 'Y', 'Z', 'True_Cluster_ID']].head())
```

## Phase 3: Network Construction
We need to pass the coordinates of each galaxy (node) to the network. We also need a linking length ($r_c$). If the distance between two galaxies is smaller than this length, they will be connected via an edge.

For this purpose, we take 3 steps: 

1. **Form our tree of nodes** using a KDTree for efficiency.

2. **Identify the linking length ($r_c$)** by finding the percolation threshold. We search the range of $1$ Mpc $\leq r_c \leq 25$ Mpc.
   * We plot the size of the Giant Connected Component (GCC) against $r_c$.
   * We use a **heuristic selection rule**: choosing the $r_c$ where the GCC growth rate (gradient) is maximized. This identifies the "phase transition" where the cosmic web structure begins to dominate the topology.

3. **Form the final network graph** using the chosen $r_c$.

```python
import networkx as nx
import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt

# build a KDTree for efficient neighbor searching
coordinates = my_galaxies[['X', 'Y', 'Z']].values
tree = KDTree(coordinates)
number_of_galaxies = len(my_galaxies)

# Search for rc across a broader range to capture the transition
rc_values = np.arange(1.0, 25.0, 0.5)
gcc_coverage_percentages = []

print("Scanning linking lengths for percolation threshold...")
for rc in rc_values:
    edges = tree.query_pairs(rc)
    temp_graph = nx.Graph()
    temp_graph.add_nodes_from(range(number_of_galaxies))
    temp_graph.add_edges_from(edges)

    if len(temp_graph.edges) > 0:
        largest_cluster = max(nx.connected_components(temp_graph), key=len)
        gcc_coverage = len(largest_cluster)
    else:
        gcc_coverage = 1

    gcc_coverage_percentage = (gcc_coverage / number_of_galaxies) * 100
    gcc_coverage_percentages.append(gcc_coverage_percentage)

# Selection Rule: Maximum Gradient (Heuristic Estimate)
gradients = np.gradient(gcc_coverage_percentages)
r_c = rc_values[np.argmax(gradients)]
print(f"Heuristic connection radius (rc): {r_c} Mpc")
```

```python
# Plotting the GCC coverage
plt.figure(figsize=(10, 6))
plt.plot(rc_values, gcc_coverage_percentages, marker='o', linestyle='-', color='#45A29E', linewidth=2)
plt.axvline(x=r_c, color='#66FCF1', linestyle='--', label=f'Percolation Threshold ({r_c} Mpc)')

plt.title("Network Phase Transition: Finding the Cosmic Web", fontsize=14, color='white')
plt.xlabel("Linking Length (Mpc)", fontsize=12, color='white')
plt.ylabel("Size of Giant Component (% of Galaxies)", fontsize=12, color='white')
plt.grid(color='#1F2833', linestyle='--', linewidth=0.5)
plt.legend()

# Dark theme styling
plt.gca().set_facecolor('#0B0C10')
plt.gcf().patch.set_facecolor('#050608')
plt.tick_params(colors='white')

plt.show()
```

## Phase 4: Community Detection and Validation

```python
from networkx.algorithms import community as nx_comm
from sklearn.metrics import normalized_mutual_info_score
#Build the Graph using the optimal linking length
pairs = tree.query_pairs(r= r_c)

# Initialize the NetworkX Graph
G = nx.Graph()
G.add_nodes_from(range(len(my_galaxies)))
G.add_edges_from(pairs)

print(f"Network built! Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")

# Run Louvain Community Detection
print("\nRunning Louvain Algorithm...")
communities = nx_comm.louvain_communities(G)
modularity = nx_comm.modularity(G, communities)

print(f"Superclusters (Communities) found: {len(communities)}")
print(f"Modularity Score (Q): {modularity:.3f}")

# Map labels back to DataFrame
predicted_labels = [0] * len(my_galaxies)
for comm_id, comm in enumerate(communities):
    for node in comm:
        predicted_labels[node] = comm_id
my_galaxies['Louvain_Community_ID'] = predicted_labels

# ==========================================
# 3. Symmetric Statistical Validation
# ==========================================
print("\nCalculating Symmetric Statistical Validation...")

# To address Issue 3: Symmetric Validation
# We filter both ground-truth and Louvain partitions to include only galaxies
# belonging to clusters of size >= 3. This removes "field" galaxies and 
# algorithmic noise (singletons) from both sides for a fair comparison.

min_size = 3
catalog_group_sizes = my_galaxies.groupby('True_Cluster_ID')['objID'].transform('count')
louvain_group_sizes = my_galaxies.groupby('Louvain_Community_ID')['objID'].transform('count')

# Mask for galaxies in significant clusters on both sides
# (Excluding ID 0 which is the catalog's explicit field galaxy label)
symmetric_mask = (catalog_group_sizes >= min_size) & \
                 (louvain_group_sizes >= min_size) & \
                 (my_galaxies['True_Cluster_ID'] != 0)

nmi_score = normalized_mutual_info_score(
    my_galaxies.loc[symmetric_mask, 'True_Cluster_ID'],
    my_galaxies.loc[symmetric_mask, 'Louvain_Community_ID']
)

print(f"Symmetric NMI (Clusters size >= {min_size}): {nmi_score:.3f}")
print(f"Galaxies included in validation: {symmetric_mask.sum()} / {len(my_galaxies)}")
```

### Phase 4.1: Results Summary

```python
import pandas as pd

summary_data = {
    "Metric": [
        "Number of Galaxies (Nodes)",
        "Redshift Range",
        "Absolute Magnitude Limit (M)",
        "Linking Length (rc)",
        "Number of Edges",
        "Number of Louvain Communities",
        "Modularity (Q)",
        "Symmetric NMI (size >= 3)"
    ],
    "Value": [
        len(my_galaxies),
        f"{z_min} - {z_max}",
        f"{M_limit:.2f}",
        f"{r_c:.2f} Mpc",
        G.number_of_edges(),
        len(communities),
        f"{modularity:.3f}",
        f"{nmi_score:.3f}"
    ]
}

summary_df = pd.DataFrame(summary_data)
summary_df
```

## Phase 5: Topological Analysis

```python
print("Initiating Phase 5: Topological Analysis...")

# 1. Calculate Degree
degrees = [G.degree(n) for n in G.nodes()]
my_galaxies['Degree'] = degrees

# 2. Isolate the Cosmic Hubs
top_hubs = my_galaxies.sort_values(by='Degree', ascending=False).head(10)

print("\n--- Top 10 Cosmic Hubs Identified ---")
print(top_hubs[['ra', 'dec', 'z', 'Degree', 'Louvain_Community_ID']])

# 3. Degree Distribution (Exploratory Analysis)
degree_counts = np.bincount(degrees)
degrees_x = np.nonzero(degree_counts)[0]
counts_y = degree_counts[degrees_x]

fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
ax.scatter(degrees_x, counts_y, marker='o', color='black', s=15, alpha=0.7)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel('Degree $k$ (Number of Connections)', fontsize=12)
ax.set_ylabel('Frequency $P(k)$', fontsize=12)
ax.set_title('Topological Signature: Exploratory Degree Distribution', fontsize=12, pad=10)

plt.tight_layout()
plt.show()
```
