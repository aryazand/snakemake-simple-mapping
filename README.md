# Fork of Snakemake-simple-mapping 

This is a fork of [snakemake-simple-mapping](https://github.com/MPUSP/snakemake-simple-mapping) that has added support for: 
* **UMIs** (in process): extract UMIs and use UMIs to deduplicate
* **trim_galore**: use trim_galore to process fastq files  

## Usage

The usage of the original workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/MPUSP/snakemake-simple-mapping).

Detailed information about input data and workflow configuration can also be found in the [`config/README.md`](config/README.md).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

_Workflow overview:_

<!-- include overview-->
<img src="resources/images/dag.png" align="center" />

## Deployment options

To run the workflow from command line, change the working directory.

```bash
cd path/to/snakemake-simple-mapping
```

Adjust options in the default config file `config/config.yml`.
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run
```

To run the workflow with test files using **conda**:

```bash
snakemake --cores 2 --sdm conda --directory .test
```
