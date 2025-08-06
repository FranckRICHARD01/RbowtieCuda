# RbowtieCuda package: an R wrapper for nvBowtie, a fork version of the [NVBIO](https://nvlabs.github.io/nvbio) library. 

The `RbowtieCuda` package provides an R interface to [nvBowtie](https://github.com/NVlabs/nvbio)^[Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.], a CUDA-based adaptation of the popular Bowtie2 short-read aligner, originally developed by Jacopo Pantaleoni and collaborators. This package allows users to build genome indexes and create alignment files in .sam or .bam formats, leveraging the computational power of modern NVIDIA GPUs for significantly faster processing speeds. nvBowtie is part of the `NVBIO` library, and the package also includes an experimental implementation of the `Wavefront Alignment (WFA)` method ^[Marco-Sola S, Moure JC, Moreto M et al. Fast gap-affine pairwise alignment using the wavefront algorithm. Bioinformatics 2021;37: 456–63.].

## Source Package Information

The `RbowtieCuda` is built on a modified version of the nvBIO v1.1.50 source code, which can be accessed at https://github.com/NVlabs/nvbio. 


## Acknowledgements

We would like to thank Ismael Galve Roperh for his assistance.

Dedicated to the memory of Lilit Tabirian


## Contributors

Original NVBIO developers:

    Jacopo Pantaleoni - jpantaleoni@nvidia.com
    Nuno Subtil - nsubtil@nvidia.com

RbowtieCuda developers:

    Samuel Simon Sanchez - samsimon@ucm.es
    Franck RICHARD - franck.richard@winstars.net

Maintainer: [Franck RICHARD](mailto:franck.richard@winstars.net)


## References

[1] Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.

[2] Marco-Sola S, Moure JC, Moreto M et al. Fast gap-affine pairwise alignment using the wavefront algorithm. Bioinformatics 2021;37: 456–63.
