These are scripts for running the Artic sequencing pipeline using the [Argo Workflows tool](https://github.com/argoproj).

The file artic.yaml contains the template for running the Artic pipeline itself. The templates in final/ were for reassembling genomes for instances where we sequenced a sample multiple times for higher coverage. The scripts/ folder contains the python scripts which run during each step of the pipeline.

The scripts are configured to run using data in our lab's private Google Cloud Storage bucket and also rely on a Docker image containing sequencing tools. We are unfortunately unable to share the Dockerfile due to restrictions by Oxford Nanopore.

