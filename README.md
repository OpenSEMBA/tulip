# Tulip
[![Windows | Build and tests](https://github.com/OpenSEMBA/tulip/actions/workflows/windows-builds-and-tests.yml/badge.svg)](https://github.com/OpenSEMBA/tulip/actions/workflows/windows-builds-and-tests.yml)

[![Ubuntu | Build and tests](https://github.com/OpenSEMBA/tulip/actions/workflows/ubuntu-builds-and-tests.yml/badge.svg)](https://github.com/OpenSEMBA/tulip/actions/workflows/ubuntu-builds-and-tests.yml)

**Tulip** (**T**ransmission line per **u**nit **l**ength and **i**n-cell **p**arameters) is a solver to obtain the _per unit length_ (PUL) $C$ and $L$ matrices which characterize electromagnetic propagation within a multiconductor tranmission line (MTL). Tulip is based on finite element methods to solver an electroestatic problem for each conductor. Some of its features are:
- Calculation of p.u.l $C$ and $L$ matrices.
- Third order isoparametric elements. 
- Support for dielectric materials.
- Open boundary conditions.
- Works on closed and open MTLs.
- Multilevel domain decomposition.
- Uses a modified [MFEM](https://mfem.org/) solver engine available [here](https://github.com/OpenSEMBA/mfem).
- Result visualization with [Paraview](https://www.paraview.org/).
- Start from `.step` CAD files with the integrated [gmsh](https://gmsh.info) workflow.

## Compile and testing
Compilation needs vcpkg with the packages stated in the ```vcpkg.json``` manifest. 

Additionally needs:
- mfem (with the version pointed by the external/mfem-geg submodule)

Configure and build presets are available. For instance, to configure in windows you can use

```shell 
    cmake 
        -DCMAKE_FIND_USE_PACKAGE_REGISTRY=FALSE  
        --preset "msbuild-vcpkg"
        -S <project folder>
        -B <build folder>
```

Once compiled, test cases can be launched from the project root folder, with

```shell
   ctest --test-dir ./<build folder> --output-on-failure
```

Most cases will store results in the `Results` folder. 
Please check the codes in `test` folder for information on the validation cases and their expected tolerances.

### Container commands

These commands require a working Docker installation and permission to create
directories beside the input file. Run `make help` to see the commands from the
terminal.

Build the reusable container images and export the runnable bundle to `dist`:

```shell
make build
```

This creates the runtime image `tulip:latest`, the test image
`tulip-test:latest`, and the standalone binary bundle in `dist/`.

Create a compressed Ubuntu 26.04 release bundle. By default it uses the
current Git tag or commit as its version; set `VERSION` to provide the release
version explicitly:

```shell
make release-ubuntu VERSION=1.2.3
```

The resulting archive is written to
`release/tulip-1.2.3-ubuntu-26.04.tar.gz`.

Run the complete test suite from the prebuilt image:

```shell
make test
```

Run an input file in the container. A new results directory is created beside
the input file and its path is printed when execution finishes. The input's
directory is copied inside the container so relative CAD and mesh references
continue to work.

```shell
make run FILE=testData/empty_coax/empty_coax.tulip.input.json
```

For paths without spaces, the input can also be passed positionally:

```shell
make run testData/empty_coax/empty_coax.tulip.input.json
```

The generated directory is named `<case>.tulip-output.XXXXXX` and contains all
Tulip output, including the results JSON and the `Results/` directory. After
`make build`, both `make test` and `make run` use the already compiled images
and do not rebuild Tulip.

If the source code, Dockerfile, dependencies, or test fixtures change, run
`make build` again to refresh the images before running tests or analyses.

## Usage example
Call `tulip` from command line as,

```shell
    tulip.exe -i <input file>
```

The input file format is [described here](docs/tulip_data_format.md).

By default, `tulip` will generate a file called `matrices.tulip.out.json` which contains the $C$ and $L$ p.u.l parameters of the MTL. Each row and column corresponds to the `N` integer in `Conductor_N`. `Conductor_0` is used as reference. 
These results have been cross-compared [here][test/DriverTest.cpp] to match with [Ansys Maxwell](https://www.ansys.com/products/electronics/ansys-maxwell). 
Comparison with [SACAMOS](https://www.sacamos.org/) does not produce satisfactory because of the different underlying analytical assumptions that are made.

This means two results for each conductor different from zero: with and without accounting for dielectrics, used to compute the p.u.l $C$ and $L$ matrices, respectively.
Below there is an example of the electric fields for the `five_wires` case visualized in Paraview with (above) and without (below) considering dielectrics.

![Electric field for the five wires case with dielectrics](docs/fig/five_wires_conductor_2_with_dielectrics_E_field.png)
![Electric field for the five wires case without dielectrics](docs/fig/five_wires_conductor_2_without_dielectrics_E_field.png)

## License and copyright
``` tulip ``` is licensed under [GPL v2](LICENSE). Its copyright belongs to the University of Granada. 

## Acknowledgements
This project is funded by the following grants:

- HECATE - Hybrid ElectriC regional Aircraft distribution TEchnologies. HE-HORIZON-JU-Clean-Aviation-2022-01. European Union.
- ESAMA - Metodos numericos avanzados para el analisis de materiales electricos y magneticos en aplicaciones aerospaciales. PID2022-137495OB-C31. Spain.
