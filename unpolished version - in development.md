# Network-Based Reconstruction of Galaxy Structure (In Development)

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

```python
import networkx as nx
import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
# build a KDTree for efficient neighbor searching
coordinates = my_galaxies[['X', 'Y', 'Z']].values
tree = KDTree(coordinates)
number_of_galaxies = len(my_galaxies)

#find the rc to determine edge connections
rc_values = np.arange(.1, 12, .2)
gcc_coverage_percentages = []

for rc in rc_values:
    # find neighbors within the connection radius
    edges = tree.query_pairs(rc)

    # create a temporary graph and add edges
    temp_graph = nx.Graph()
    temp_graph.add_nodes_from(range(number_of_galaxies))
    temp_graph.add_edges_from(edges)

    # find the largest connected component
    if len(temp_graph.edges) > 0:
        largest_cluster = max(nx.connected_components(temp_graph), key=len)
        gcc_coverage = len(largest_cluster)
    else:
        gcc_coverage = 1  # if no edges, the largest cluster is just one galaxy

    # calculate the percentage of galaxies in the largest cluster
    gcc_coverage_percentage = (gcc_coverage / number_of_galaxies) * 100
    gcc_coverage_percentages.append(gcc_coverage_percentage)

# finding the optimal rc using the maximum gradient method
gradients = np.gradient(gcc_coverage_percentages)
r_c = rc_values[np.argmax(gradients)]
print("Optimal connection radius (rc) for percolation threshold:", r_c)
```

```python
# Visualizing the percolation threshold
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

# Run Louvain Algorithm
print("\nRunning Louvain Algorithm...")
communities = nx_comm.louvain_communities(G)
modularity = nx_comm.modularity(G, communities)

print(f"Superclusters (Communities) found: {len(communities)}")
print(f"Modularity Score (Q): {modularity:.3f}")

# Statistical Validation
predicted_labels = [0] * len(my_galaxies)
for comm_id, comm in enumerate(communities):
    for node in comm:
        predicted_labels[node] = comm_id
my_galaxies['Louvain_Community_ID'] = predicted_labels

true_labels = my_galaxies['True_Cluster_ID'] != 0
nmi_score = normalized_mutual_info_score(my_galaxies.loc[true_labels, 'True_Cluster_ID'],
                                          my_galaxies.loc[true_labels, 'Louvain_Community_ID'])

print(f"Normalized Mutual Information (NMI): {nmi_score:.3f}")
```

## Phase 5: Visualizing Filaments and Nodes

```python
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

print("Generating publication-ready filament plot...")

fig, ax = plt.subplots(figsize=(14, 7), dpi=150)
plt.rcParams.update({'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': 1.2})

# Plot Background
ax.scatter(my_galaxies['X'], my_galaxies['Y'], s=0.5, color='black', alpha=0.5, zorder=1)

# Plot Massive Nodes
top_cluster_ids = my_galaxies['Louvain_Community_ID'].value_counts().head(75).index
cluster_centers_x = []
cluster_centers_y = []

for cid in top_cluster_ids:
    cluster_data = my_galaxies[my_galaxies['Louvain_Community_ID'] == cid]
    cluster_centers_x.append(cluster_data['X'].mean())
    cluster_centers_y.append(cluster_data['Y'].mean())

ax.scatter(cluster_centers_x, cluster_centers_y, s=70, color='red', zorder=3)

# Draw Filaments
edges_list = list(G.edges())
np.random.seed(42)
sampled_indices = np.random.choice(len(edges_list), 15000, replace=False)

lines = []
for idx in sampled_indices:
    u, v = edges_list[idx]
    x0, y0 = my_galaxies.iloc[u]['X'], my_galaxies.iloc[u]['Y']
    x1, y1 = my_galaxies.iloc[v]['X'], my_galaxies.iloc[v]['Y']
    lines.append([(x0, y0), (x1, y1)])

lc = LineCollection(lines, colors='#3b82f6', linewidths=0.8, alpha=0.6, zorder=2)
ax.add_collection(lc)

ax.set_xlabel('X [Mpc]', fontsize=14)
ax.set_ylabel('Y [Mpc]', fontsize=14)
ax.set_title('2D Projected Cosmic Web: Filaments and Nodes', fontsize=14, pad=10)
ax.set_aspect('equal')

plt.tight_layout()
plt.show()
```

## Phase 6: Richness Distribution Analysis

In this phase, we compare the **richness** (the number of galaxies per cluster) of our detected Louvain communities against the ground truth from the Tempel et al. (2017) catalog.

### Why this matters:
1. **Physical Scale Validation**: It verifies if our network-based clusters match the expected size distribution of real astronomical groups.
2. **Linking Length Audit**: It helps determine if our chosen linking length ($r_c$) is physically appropriate. If our clusters are systematically much larger or smaller than the catalog, it suggests we may be over-linking or under-linking the data.
3. **Scaling Behavior**: Cosmic structures typically follow a power-law distribution in richness. We will visualize this on a log-log plot to see if our network captures this fundamental scaling property of the Cosmic Web.

```python
import seaborn as sns

print("Calculating richness distributions...")

# 1. Calculate richness
louvain_richness = my_galaxies['Louvain_Community_ID'].value_counts().values
catalog_richness = my_galaxies[my_galaxies['True_Cluster_ID'] != 0]['True_Cluster_ID'].value_counts().values

# 2. Plotting
plt.figure(figsize=(10, 6), dpi=150)
sns.histplot(catalog_richness, label='Tempel et al. (2017) Catalog', element="step", fill=False, color='#F5BDE6', log_scale=(True, True), linewidth=2)
sns.histplot(louvain_richness, label='Louvain Communities (Predicted)', element="step", fill=False, color='#8AADF4', log_scale=(True, True), linewidth=2)

plt.title("Richness Distribution: Catalog vs. Network Communities", fontsize=14, color='white')
plt.xlabel("Richness (Galaxies per Cluster)", fontsize=12, color='white')
plt.ylabel("Frequency (Count)", fontsize=12, color='white')
plt.legend()

plt.gca().set_facecolor('#0B0C10')
plt.gcf().patch.set_facecolor('#050608')
plt.tick_params(colors='white')
plt.grid(color='#1F2833', linestyle='--', linewidth=0.5, alpha=0.3)

plt.show()

# Summary Statistics
print("\n--- Richness Summary Statistics ---")
print(f"Average Catalog Richness: {catalog_richness.mean():.2f}")
print(f"Average Louvain Richness: {louvain_richness.mean():.2f}")
print(f"Max Catalog Richness: {catalog_richness.max()}")
print(f"Max Louvain Richness: {louvain_richness.max()}")
```
