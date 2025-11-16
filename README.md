<h1 align="center">xCell TME Enrichment – Python rpy2 Wrapper</h1>

<p align="center">
  A Python 3.12 command-line wrapper around the R-based <strong>xCell</strong> pipeline for
  tumor microenvironment (TME) <strong>cell-type enrichment</strong> profiling.
</p>

<hr />

<h2>Overview</h2>

<p>
  This repository provides a reproducible Python interface to the
  <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1" target="_blank" rel="noreferrer">
    xCell method (Aran et&nbsp;al., Genome Biology 2017)
  </a>,
  which estimates enrichment of 64 immune and stromal cell types from bulk gene expression data.
</p>

<p>
  The script <code>xcell_tme_deconvolution_runner.py</code> is a thin but robust wrapper around an
  R function (<code>run_xcell_pipeline()</code>) implemented via <code>rpy2</code>. It:
</p>

<ul>
  <li>Runs the full xCell enrichment and spillover compensation pipeline.</li>
  <li>Supports multiple expression file formats (<code>.csv</code>, <code>.tsv</code>, <code>.txt</code>, <code>.xlsx</code>).</li>
  <li>Automatically detects whether genes are in rows or columns and re-orients as needed.</li>
  <li>Handles Ensembl-to-HGNC symbol mapping and duplicates collapsing.</li>
  <li>Applies data-type-aware normalization (edgeR TMM, log transforms, or pass-through).</li>
  <li>Produces xCell scores, composite TME indices, heatmaps, PCA/UMAP plots, and stacked compositions.</li>
  <li>Writes biological context and automated interpretations for each cell type.</li>
</ul>

<hr />

<h2>Repository Layout</h2>

<p>Typical layout for this project:</p>

<ul>
  <li><code>xcell_tme_deconvolution_runner.py</code> – main Python CLI wrapper.</li>
  <li><code>requirements.txt</code> – Python dependencies (including <code>rpy2</code>).</li>
  <li><code>Xcell_Paper_ref.pdf</code> – reference PDF of the original xCell paper.</li>
  <li><code>Teset_Run/</code> – example folder with test data and outputs.</li>
  <li><code>.venv/</code> – optional local virtual environment (recommended).</li>
</ul>

<p>
  You can extend this with more structured folders (e.g. <code>xcell_enrichment/</code>, <code>xcell_tests/</code>,
  <code>xcell_docs/</code>) as the project grows.
</p>

<hr />

<h2>Requirements</h2>

<h3>Python</h3>
<ul>
  <li><strong>Python 3.12</strong> (tested)</li>
  <li>Recommended OS: Windows, Linux, or macOS with a working R installation.</li>
  <li>Key Python package: <code>rpy2</code> (installed via <code>requirements.txt</code>).</li>
</ul>

<h3>R</h3>
<ul>
  <li><strong>R &gt;= 4.1</strong> (recommended)</li>
  <li>The script can optionally auto-install R dependencies when run with <code>--do-install</code>:</li>
</ul>

<pre><code>CRAN:
  data.table, readr, stringr, edgeR, pheatmap, ggplot2, RColorBrewer, viridis,
  dplyr, tidyr, tibble, knitr, kableExtra, rmarkdown, uwot, readxl,
  cowplot, gridExtra

Bioconductor:
  xCell, AnnotationDbi, org.Hs.eg.db
</code></pre>

<p>
  If you prefer to manage R packages yourself (e.g. on HPC), you can install these manually
  and run the pipeline with <code>--no-install</code>.
</p>

<hr />

<h2>Python Virtual Environment (recommended)</h2>

<p>Create and use a dedicated virtual environment for this project.</p>

<h3>Windows (PowerShell)</h3>

<pre><code>cd &lt;path-to-your-repo&gt;

python -m venv .venv
.\.venv\Scripts\Activate.ps1

pip install --upgrade pip
pip install -r requirements.txt
</code></pre>

<h3>Linux / macOS</h3>

<pre><code>cd &lt;path-to-your-repo&gt;

python3 -m venv .venv
source .venv/bin/activate

pip install --upgrade pip
pip install -r requirements.txt
</code></pre>

<hr />

<h2>Input Files</h2>

<h3>1. Expression matrix (<code>--expr-file</code>)</h3>

<p>The expression table can be:</p>

<ul>
  <li>Genes in rows, samples in columns, first column = gene ID.</li>
  <li>Or samples in rows, genes in columns, first column = sample ID.</li>
</ul>

<p>The pipeline will:</p>

<ul>
  <li>Detect the orientation and transpose if necessary.</li>
  <li>Convert numeric values appropriately and handle non-numeric entries gracefully.</li>
  <li>Optionally map Ensembl IDs (e.g. <code>ENSG...</code>) to HGNC symbols via <code>org.Hs.eg.db</code>.</li>
  <li>Collapse duplicated gene symbols by summation.</li>
</ul>

<h3>2. Metadata table (<code>--meta-file</code>)</h3>

<p>
  The metadata table must contain at least the following columns (case-insensitive,
  the script normalizes header names):
</p>

<ul>
  <li><code>sample_id</code> – raw sample identifier matching the expression matrix columns/rows.</li>
  <li><code>condition</code> – grouping variable (e.g. <code>AKI</code>, <code>Control</code>, etc.).</li>
</ul>

<p>
  Internally, the pipeline standardizes sample IDs (lowercase, strips extensions, cleans symbols)
  and aligns metadata with expression samples. Samples present in metadata but not in the expression
  matrix are dropped with a note in the logs.
</p>

<hr />

<h2>Running the Pipeline</h2>

<h3>Basic command</h3>

<p>
  Once your virtual environment is activated and dependencies are installed, you can run:
</p>

<pre><code>python xcell_tme_deconvolution_runner.py \
  --data-dir "D:\AyassBio_Workspace_Downloads\xcell-deonv\xcell_py_pkg\Teset_Run" \
  --expr-file "GSE139061_Eadon_processed_QN_101419.csv" \
  --meta-file "sample_metadata_patientAKI.tsv" \
  --out-dir "xcell_output_V31_BIO" \
  --no-install
</code></pre>

<p>
  Notes:
</p>

<ul>
  <li><code>--data-dir</code> is used as a base path to resolve relative <code>--expr-file</code>, <code>--meta-file</code>, and <code>--out-dir</code>.</li>
  <li>You can provide absolute paths for any of these arguments; they will be respected.</li>
  <li><code>--no-install</code> skips automatic R package installation (recommended if your R environment is already set up).</li>
</ul>

<h3>Key CLI arguments</h3>

<ul>
  <li><code>--data-dir</code> (optional): base directory for resolving relative paths.</li>
  <li><code>--expr-file</code> (required): expression matrix (csv/tsv/txt/xlsx).</li>
  <li><code>--meta-file</code> (required): metadata table with <code>sample_id</code> and <code>condition</code>.</li>
  <li><code>--out-dir</code> (required): output directory for all results.</li>
  <li><code>--do-install</code> / <code>--no-install</code>: whether to install/ensure R dependencies from within the pipeline.</li>
  <li><code>--thresh-med-high</code>: median threshold for calling enrichment “High” (default 0.20).</li>
  <li><code>--thresh-med-mod</code>: median threshold for calling enrichment “Moderate” (default 0.10).</li>
  <li><code>--pres-fdr-alpha</code>: FDR cut-off for presence significance (default 0.10).</li>
  <li><code>--pres-min-frac</code>: minimum fraction of significant samples to call “Frequent” (default 0.35).</li>
  <li><code>--umap-min-samples</code>: minimum number of samples required to compute UMAP (default 20).</li>
  <li><code>--max-cards</code>: maximum number of “sample cards” (not currently used in plotting; accepted for future extension).</li>
  <li><code>--top-n-card</code>: top-N cell types per sample card (currently used for limiting boxplot subsets).</li>
</ul>

<hr />

<h2>Outputs</h2>

<p>
  All outputs are written under <code>--out-dir</code>. The pipeline creates a structured set of subfolders:
</p>

<ul>
  <li><code>xCell_Scores/</code>
    <ul>
      <li><code>xcell_raw_scores.csv</code> – raw enrichment scores from xCell.</li>
      <li><code>xcell_transformed_scores.csv</code> – scores after curve fitting / transformation.</li>
      <li><code>xcell_final_scores.csv</code> – spillover-compensated scores (primary matrix to use).</li>
      <li><code>Cohort_TME_Matrix.csv</code> – per-sample summary table containing:
        <ul>
          <li><code>sample_std</code> – standardized sample ID.</li>
          <li><code>condition</code>.</li>
          <li><code>ImmuneScore_xCell</code>.</li>
          <li><code>StromaScore_xCell</code>.</li>
          <li><code>MicroenvironmentScore_xCell</code>.</li>
          <li><code>TME_Inflammation_Index</code> – “hot vs. cold” TME proxy.</li>
        </ul>
      </li>
      <li><code>xcell_cohort_celltype_enrichment_summary.csv</code> – per-cell-type median, presence fraction, and narrative interpretation.</li>
      <li><code>High_Variability_CellTypes.csv</code> – cell types with many significant samples by FDR.</li>
      <li><code>normalized_expression_&lt;method&gt;.csv</code> – normalized expression matrix (edgeR TMM logCPM, log2, or as-is).</li>
    </ul>
  </li>

  <li><code>Cell-Type Enrichment Profile (Per Sample)/</code>
    <ul>
      <li><code>TME_Enrichment_Heatmap_AllCellTypes.png</code> – heatmap of all cell types across samples.</li>
      <li><code>TopVariability_CellTypes_Heatmap.png</code> – heatmap of the most variable cell types.</li>
      <li><code>xCell_Score_PCA_by_Condition.png</code> – PCA of xCell scores colored by condition.</li>
      <li><code>CellType_Correlation_Spearman.png</code> – cell-type correlation matrix.</li>
      <li><code>UMAP_xCell_EnrichmentSpace.png</code> – UMAP of xCell enrichment space (if enough samples).</li>
      <li><code>Composite_Indices_by_Condition.png</code> – ImmuneScore, StromaScore, MicroenvironmentScore, and TME inflammation index by condition.</li>
      <li><code>Stacked_Composition_AllCellTypes.png</code> – stacked composition plot of enrichment per sample.</li>
      <li><code>Boxplots_TopK/</code> – boxplots of top-variability cell types stratified by condition.</li>
      <li><code>Xcell_Analytical_Figures/</code> – reserved for additional figures/cards (future extension).</li>
    </ul>
  </li>

  <li><code>ImmuneProfile_Context/</code>
    <ul>
      <li><code>xcell_bio_context.csv</code> – high-level biological descriptions for each xCell cell type and composite score.</li>
    </ul>
  </li>

  <li><code>sessionInfo.txt</code> – full R session information for reproducibility (R version, packages, etc.).</li>
</ul>

<hr />

<h2>Reproducibility & Design Choices</h2>

<ul>
  <li>All key steps are logged via intermediate CSV files and plots.</li>
  <li>R session info is saved to <code>sessionInfo.txt</code> to capture exact package versions.</li>
  <li>Orientation detection for expression matrices is automatic and debug CSVs are written if alignment fails.</li>
  <li>All file paths used by R are passed explicitly from the Python CLI; no hidden defaults are used for I/O.</li>
</ul>

<hr />

<h2>Troubleshooting</h2>

<h3><code>ImportError: No module named 'rpy2'</code></h3>
<ul>
  <li>Ensure your virtual environment is activated.</li>
  <li>Re-run <code>pip install -r requirements.txt</code>.</li>
</ul>

<h3>R package installation issues</h3>
<ul>
  <li>If automatic installation (<code>--do-install</code>) fails, install the listed CRAN/Bioconductor packages manually in R, then run with <code>--no-install</code>.</li>
  <li>Check that your R installation is on the system PATH and compatible with your Python/rpy2 build.</li>
</ul>

<h3>No overlapping sample IDs between metadata and expression</h3>
<ul>
  <li>Verify that <code>sample_id</code> in the metadata matches the expression sample IDs.</li>
  <li>Remember that the script standardizes IDs (lowercases, strips extensions, and cleans symbols).</li>
  <li>Inspect generated debug files <code>_debug_plotmat_cols.csv</code> and <code>_debug_requested_samples.csv</code> in the output directory.</li>
</ul>

<hr />

<h2>Citation</h2>

<p>
  If you use this wrapper in your work, please cite both the original xCell paper and this repository:
</p>

<ul>
  <li>
    <strong>xCell method</strong>:<br />
    Aran D, Hu Z, Butte AJ. <em>xCell: digitally portraying the tissue cellular heterogeneity landscape.</em>
    Genome Biology. 2017;18:220. doi:10.1186/s13059-017-1349-1.
  </li>
  <li>
    <strong>Python wrapper</strong>:<br />
    Please reference this GitHub repository (URL and commit) in the Methods section as the implementation used for running xCell.
  </li>
</ul>

<hr />
