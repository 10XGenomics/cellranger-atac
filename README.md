# Building Cell Ranger ATAC 1.2.0
## Build dependencies

This package has been built internally using:
- Python 2.7.13
- clang 7.0 (gcc/g++ should also work)
- go 1.12

## Setting up the environment and build dependencies

- Environment dependencies can be found in the Ranger ATAC 1.2.0 package (https://support.10xgenomics.com/developers/software/downloads/latest)
  - The Ranger package includes a build of the Martian platform (v3.2.4), which is open source. For more information, go to (http://martian-lang.org).
  - The Ranger package includes a build of miniconda, which contains required python executables and packages.

### Example setup of build dependencies on Ubuntu 20.04
```
# use gcc/g++ in this case to support GLIBCXX_3.4.25 in pre-built ranger package
sudo apt-get install make gcc g++ pkg-config

# Add golang to path
wget https://golang.org/dl/go1.12.17.linux-amd64.tar.gz
tar xf go1.12.17.linux-amd64.tar.gz
export PATH=/path/to/go/bin:$PATH
export GOROOT=/path/to/go

# Source the Ranger package for python, the Martian Environment and binary dependencies
source /path/to/ranger-1.2.0-ATAC/sourceme.bash
```

## Build command
```
make
```

# Running Cell Ranger ATAC
## Setting up the Runtime dependencies for Cell Ranger ATAC
```
source /path/to/open-source/cellranger-atac/sourceme.bash
```

## Note about Loupe
The binaries required to generate Loupe Cell Browser (.cloupe) are not included in this repository or in the binary dependencies package Ranger. The necessary binaries, including crconverter, bamtofastq, etc. can be obtained from an existing binary version of Cell Ranger ATAC by running:
```
cp /path/to/cellranger-atac/cellranger-atac-cs/1.2.0/lib/bin/crconverter /path/to/open-source/cellranger-atac/lib/bin/
```

# Support
We do not provide support for building and running this code.

The officially supported release binaries are available at: (https://support.10xgenomics.com/single-cell-atac/software/downloads/latest)

