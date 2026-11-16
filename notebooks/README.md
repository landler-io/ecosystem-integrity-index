# Notebooks

Example notebooks for the EII client.

| Notebook | Content |
|---|---|
| `eii_client_tutorial.ipynb` | Point values, polygon statistics, zonal statistics, rasters, Earth Engine layers, change between reference years |
| `local_modulation_tutorial.ipynb` | Natural Capital modulation of EII with default, pre-computed or custom KPIs |

Both set `REFERENCE_YEAR = 2024` at the top. Without a `reference_year`, the client uses the
latest year (2025).

## Running the notebooks

1. Install the package with the client extras and the plotting dependencies of `quicklook`:

   ```bash
   pip install "ecosystem-integrity-index[client] @ git+https://github.com/landler-io/ecosystem-integrity-index"
   pip install jupyter matplotlib contextily dask
   ```

2. Authenticate with Earth Engine: `earthengine authenticate`.
3. Replace `your-gee-project-id` in the setup cell with your Cloud project and run the notebook
   in Jupyter or VS Code.
