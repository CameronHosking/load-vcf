# load-vcf
Load a vcf file and quickly separate it by sample. Mostly a helper file for other projects.

#### Sample binary vcf format
This format supports: 
* Up to length 255 variant IDs, and allelels. 
* Haploid and Diploid and phased or unphased. 
* Missing values in the vcf (.) which are treated as referance values.

The resulting .bvcf file is a concatination of all variants in a sample.
Each sample is encoded as following

Flag bit set |  Meaning
------------- | -------------
128  |  set for Diploid unset for haploid
64 | set for homozygous unset for heterozygous
32 | set for phased unset for unphased
16 | set if first ploid is the referance allele
8 | set if the second ploid is the referance allele

Location |  Content
------------- | -------------
Bytes 1-4  | Variant position in genome (little endian 32 bit uint)
Byte 5 | Flag Byte (see above)
Byte 6 | Length of referance sequence (8 bit uint)
Byte 7 | Length of variant ID (8 bit uint)
Byte 8 | Length of first allele (0 if referance) (8 bit uint)
Byte 9 | Length of second allele (0 if referance or homozygous) (8 bit uint)

Then follows a string of length Byte 7 + Byte 8 + Byte 9 which is the Variant ID followed by the first allele then the second allele. These are not null terminated strings.
