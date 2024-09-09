# 🧬 Bioinformatics on Flyte
---

This repo contains tasks, workflows, image definitions, and datatypes used to standardize the orchestration of common bioinformatics tasks using Flyte.

# Quickstart

1. Spin up the local development [sandbox](https://docs.flyte.org/en/latest/user_guide/environment_setup.html#create-a-local-demo-flyte-cluster)
2. `pip install union`
3. `git clone https://github.com/unionai-oss/unionbio.git && cd unionbio`
4. `union --config ~/.flyte/config-sandbox.yaml run --remote src/unionbio/workflows/simple_variant_calling.py calling_wf`

    Navigate to the local URL that's produced and watch your workflow run!

# Highlights

## 🐳 Container Images
✨**Adding custom dependencies alongside Flytekit**

[ImageSpecs](https://docs.flyte.org/projects/cookbook/en/latest/auto_examples/customizing_dependencies/image_spec.html#image-spec-example) contained in the [images](images.py) module build a standard set of OCI-compliant container images for use throughout the different workflows. They can be built with entrypoints present in [pyproject.toml](pyproject.toml).

## 📈 Datatypes
✨**Leverage dataclasses to keep things organized**

Using [dataclasses](src/unionbio/datatypes/) to define your samples provides a clean and extensible data structure to keep your workflows tidy. Instead of writing to directories and keeping track of things manually on the commandline, these dataclasses will capture relevant metadata about your samples and let you know where to find them in object storage.

## 🔍 Quality Control and Pre-processing

### FastQC
✨**Run arbitrary shell commands**

FastQC is a very common tool written in Java for gathering QC metrics about raw reads. It doesn't have any python bindings, but luckily Flyte lets us run arbitrary [ShellTasks](src/unionbio/tasks/fastqc.py) with a clean way of passing in inputs and receiving outputs. Just define a script for what you need to do and ShellTask will handle the rest.

### Automatic QC checkpointing
✨**Decide wether to continue workflow execution based on QC metrics via conditionals**

FastQC generates a summary file with a simple PASS / WARN / FAIL call across a number of different metrics. We can use [conditionals](https://docs.flyte.org/projects/cookbook/en/latest/auto_examples/advanced_composition/conditions.html) in our workflow to check for any FAIL lines in the summary and automatically halt execution. This can surface an early failure without wasting valuable compute or anyone's time doing manual review.

### FastP
✨**Specify resources and parallelize via map task**

FastP is another common pre-processing tool for filtering out bad reads, trimming, and adapter removal. It can be a bit more memory hungry than Flyte's defaults are set to; luckily we can use [Resources](src/unionbio/tasks/fastp.py) to bump that up and allow it to run efficiently. Additionally, we can make use of a [map task](https://docs.flyte.org/projects/flytekit/en/latest/_modules/flytekit/core/array_node_map_task.html) in our [workflow](src/unionbio/workflows/compare_aligners.py) to parallelize the processing of fastp across all our samples.

## 👩‍🔬 Human-in-the-Loop Approval
✨**Pause processing while waiting for human input**

As a final check before moving onto the alignment, we can define an explicit approval right in the workflow. Aggregating reports of all processing done up to this point, and visualizing it via Decks (more on that later), a researcher is able to quickly get a high level view of the work done so far and approve the analysis for further processing.

## 📏 Alignment

### Generate indices
✨**Leverage caching to save time on successive runs**

Index generation can be a very compute intensive step. Luckily, we can take advantage of Flyte's native caching when building that index for [bowtie](src/unionbio/tasks/bowtie2.py) and [hisat](src/unionbio/tasks/hisat2.py). We've also defined a `cache_version` in the [config](src/unionbio/config.py) that relies on a hash of the reference location in the object store. This means that changing the reference will invalidate the cache and trigger a rebuild, while allowing you to go back to your old reference with impunity.

### Bowtie2 vs Hisat2
✨**Compare aligners across an arbitrary number of inputs via dynamic workflows**

When prototyping a new pipeline, it's usually a good idea to evaluate a few different tools to see how they perform with respect to runtime and resource requirements. This is easy with a [dynamic](https://docs.flyte.org/en/latest/user_guide/advanced_composition/dynamic_workflows.html#dynamic-workflow) workflow, which allows us to pass in an arbitrary number of inputs to be used with whatever tasks we want. In the main [workflow](src/unionbio/workflows/compare_aligners.py) you'll pass a list of filtered samples to each tool and be able to capture run statistics in the Alignment dataclass as well as visualize their runtimes in the Flyte console.

## 📋 Reporting
✨**Visualize performance via Decks**

We use [MultiQC](src/unionbio/tasks/multiqc.py), an excellent multi-modal visualization tool for reporting. After gathering all relative metrics from a workflow, we're able to render that report via [Decks](https://docs.flyte.org/projects/cookbook/en/latest/auto_examples/development_lifecycle/decks.html), giving us rich run statistics without ever leaving the Flyte console!