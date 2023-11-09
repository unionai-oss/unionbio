# Variant Discovery on Flyte

This bioinformatics workflow will take you through a variant discovery pipeline, from raw reads to actionable insights, all orchestrated with Flyte. The workflow is broken up into sensible steps, each highlighting some of the more advanced Flyte features.

## QC and Alignment

The first section will collect QC metrics, preprocess our samples, and align them to a reference genome using a couple different aligners for comparison sake. 

### Building the Base Image
✨**Adding custom dependencies alongside Flytekit**

### Logically group samples
✨**Leverage dataclasses to keep things organized**

### FastQC
✨**Run arbitrary shell commands**

### FastP
✨**Specify resources and parallelize via map task**

### Approve QC report
✨**Require a human-in-the-loop before proceeding**

### Generate indexes
✨**Leverage caching to save time on successive runs**

### Bowtie2 vs Hisat2
✨**Compare aligners across an arbitrary number of inputs via dynamic workflows**

### MultiQC
✨**Add a dependency in-line and visualize performance via Decks**