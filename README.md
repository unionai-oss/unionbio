# 🧬 Variant Discovery on Flyte
---

This bioinformatics workflow will take you through a variant discovery pipeline, from raw reads to actionable insights, all orchestrated with Flyte. The workflow is broken up into sensible steps, each highlighting some of the more advanced Flyte features.

## 🔍📏 QC and Alignment

The first section will collect QC metrics, preprocess our samples, and align them to a reference genome using a couple different aligners for comparison sake. 

### 1. Building the Base Image
✨**Adding custom dependencies alongside Flytekit**

You can build out the default Dockerfile that ships with `pyflyte init` to add custom dependencies. The [Dockerfile](Dockerfile) here incorporates a number of common bioinformatics tools to create a base image for the project. Any additional python dependencies can also be specified in the [requirements.txt](requirements.txt) file, which gets pulled into and installed in the final image. You'll also be able to use [fast registration](https://docs.flyte.org/projects/cookbook/en/latest/getting_started/package_register.html#fast-registration) as you iterate on your code to avoid having to rebuild and push an image every time.

### 2. Logically group samples
✨**Leverage dataclasses to keep things organized**

Using [dataclasses](workflows/sample_types.py) to define your samples provides a clean and extensible data structure to keep your workflows tidy. Instead of writing to directories and keeping track of things manually on the commandline, these dataclasses will keep relevant metadata about your samples and let you know where to find them in object storage.

### 3. FastQC
✨**Run arbitrary shell commands**

FastQC is a very common tool for gathering QC metrics about raw reads. It's a java tool and included in our Dockerfile. It doesn't have any python bindings, but luckily Flyte lets us run arbitrary [ShellTasks](workflows/fastqc.py) with a clean way of passing in inputs and receiving outputs. Just define a script for what you need to do and ShellTask will handle the rest.

### 4. Approve QC report
✨**Require a human-in-the-loop before proceeding**

TBD

### 5. FastP
✨**Specify resources and parallelize via map task**

FastP is another common pre-processing tool for filtering out bad reads, trimming, and adapter removal. It can be a bit more memory hungry than Flyte's defaults are set to; luckily we can use [Resources](workflows/fastp.py) to bump that up and allow it to run efficiently. Additionally, we can make use of a [map task](https://docs.flyte.org/projects/flytekit/en/latest/_modules/flytekit/core/array_node_map_task.html) in our [workflow](workflows/alignment.py) to parallelize the processing of fastp across all our samples. 

### 6. Generate indexes
✨**Leverage caching to save time on successive runs**

Index generation can be a very compute intensive step. Luckily, we can take advantage of Flyte's native caching when building that index for [bowtie](workflows/bowtie2.py) and [hisat](workflows/hisat2.py). We've also defined a `cache_version` in the [config](workflows/config.py) that relies on a hash of the reference location in the object store. This means that changing the reference will invalidate the cache and trigger a rebuild, while allowing you to go back to your old reference with impunity.

### 7. Bowtie2 vs Hisat2
✨**Compare aligners across an arbitrary number of inputs via dynamic workflows**

When prototyping a new pipeline, it's usually a good idea to evaluate a few different tools to see how they perform with respect to runtime and resource requirements. This is easy with a `dynamic` workflow, which allows us to pass in an arbitrary number of inputs to be used with whatever tasks we want. In the main [workflow](workflows/alignment.py) you'll pass a list of filtered samples to each tool and be able to capture run statistics in the SamFile dataclass as well as visualize their runtimes in the Flyte console.

### 8. MultiQC
✨**Add a dependency in-line and visualize performance via Decks**

Finally, it's easy to tack on a dependency that's not included in your Dockerfile via an [ImageSpec](https://docs.flyte.org/projects/cookbook/en/latest/auto_examples/customizing_dependencies/image_spec.html#image-spec-example) definition inline. We do this for [MultiQC](workflows/multiqc.py), an excellent multi-modal visualization tool. After gathering all relative metrics from the workflow, we're able to render that report via [Decks](https://docs.flyte.org/projects/cookbook/en/latest/auto_examples/development_lifecycle/decks.html), giving us rich run statistics without ever leaving the Flyte console!