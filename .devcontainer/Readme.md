# ESPResSo

## Getting started with Codespaces

### Starting the IDE

* go to your ESPResSo fork
* follow the GitHub documentation on [Creating a codespace for a repository](https://docs.github.com/en/codespaces/developing-in-a-codespace/creating-a-codespace-for-a-repository#creating-a-codespace-for-a-repository)

### Building ESPResSo

General workflow:

* load MPI module
* build ESPResSo
* create Python environment

Command lines:

```sh
# load MPI module
. /etc/profile.d/modules.sh
module load mpi
# build ESPResSo
mkdir build
cd build
cmake .. -D ESPRESSO_BUILD_WITH_FFTW=ON -D ESPRESSO_BUILD_WITH_WALBERLA=ON
make -j$(nproc)
# create Python environment
python -m venv espresso
. espresso/bin/activate
realpath src/python > "$(python -c 'import sysconfig;print(sysconfig.get_path("platlib"))')/espresso.pth"
```

### Running interactive notebooks

Two options:

* run tutorials from `doc/tutorials` in the top-level directory
* run `make -j$(nproc) tutorials` to generate notebooks with hidden answers and run tutorials from `doc/tutorials` in the build directory

When the IDE asks for the interpreter, install the Python + Jupyter extensions and select the existing venv environment.
