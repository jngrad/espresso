# EasyBuild workflow

This folder contains the EasyBlock and EasyConfig required
to build the development version of ESPResSo.

```sh
# create local environment
module load EasyBuild/4.9.4
mkdir -p "${HOME}/easybuild/cache"
cd "${HOME}/easybuild"

# build toolchain
eb Automake-1.16.5-GCCcore-11.3.0.eb --robot --prefix="${HOME}/easybuild/cache"

# download dependencies
wget -LO espresso-237c376.tar.gz  "https://github.com/jngrad/espresso/archive/237c376.tar.gz"
wget -LO walberla-f36fa0a6.tar.gz "https://i10git.cs.fau.de/api/v4/projects/walberla%2Fwalberla/repository/archive?sha=f36fa0a6"
wget -LO h5xx-0.9.1.tar.gz        "https://github.com/h5md/h5xx/archive/0.9.1.tar.gz"
mv espresso-237c376.tar.gz  cache/sources/e/ESPResSo/
mv walberla-f36fa0a6.tar.gz cache/sources/e/ESPResSo/
mv h5xx-0.9.1.tar.gz        cache/sources/e/ESPResSo/

# build development version of ESPResSo
eb --force --include-easyblocks ./espresso.py ./ESPResSo-4.3-foss-2023b.eb \
   --robot --prefix="${HOME}/easybuild/cache"

# use ESPResSo
module load ESPResSo/237c376-foss-2023b
```

