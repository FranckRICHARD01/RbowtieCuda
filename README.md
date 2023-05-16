# RbowtieCuda
Bioconductor package: an R wrapper for mvBowtie, a fork version of [NVBIO](https://nvlabs.github.io/nvbio). 

The `RbowtieCuda` package provides an R interface to the [nvBowtie](https://github.com/NVlabs/nvbio)[1] short read aligner adaptation of bowtie2 unning under Cuda proposed by Jacopo Pantaleoni and al. The `RbowtiCuda` package allows users to build indexes and to create alignment files (.sam or .bam) but with a higher processing speed thanks to the computing power of modern video cards. nvBowtie is a part of the NVBIO library.

## BowtieCuda Source Package Information

The `RbowtieCuda` package uses a modificated version of the nvBio v1.1.50 source code which was obtained from https://github.com/NVlabs/nvbio. 


## Acknowledgement

We would like to thank Ismael Galve Roperh for his assistance.


## Credits

The main contributors of the original NVBIO are:

    Jacopo Pantaleoni - jpantaleoni@nvidia.com
    Nuno Subtil - nsubtil@nvidia.com
    Samuel Simon Sanchez - samsimon@ucm.es

The maintainer of the RbowtieCuda package is [Franck RICHARD](mailto:franck.richard@winstars.net)


## References

[1] Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359.
