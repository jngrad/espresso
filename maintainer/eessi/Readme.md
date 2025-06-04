# Package ESPResSo on EESSI

The following instructions package the development version of ESPResSo using EasyBuild.

```sh
export BUILDDIR="${HOME}/eessi-dev"
cd "${BUILDDIR}"
module load EasyBuild/4.9.4
mkdir -p cache

wget https://archives.boost.io/release/1.83.0/source/boost_1_83_0.tar.gz
wget https://ftp.gnu.org/gnu/groff/groff-1.23.0.tar.gz
wget https://archive.mesa3d.org/mesa-23.1.9.tar.xz
mv boost_1_83_0.tar.gz cache/sources/b/Boost.MPI/boost_1_83_0.tar.gz
mv groff-1.23.0.tar.gz cache/sources/g/groff/groff-1.23.0.tar.gz
mv mesa-23.1.9.tar.xz cache/sources/m/Mesa/mesa-23.1.9.tar.xz

mkdir -p e/ESPResSo
wget -LO cache/sources/e/ESPResSo/espresso-4b9ccc4.tar.gz   https://github.com/espressomd/espresso/archive/4b9ccc4.tar.gz
wget -LO cache/sources/e/ESPResSo/walberla-59c9b8b1.tar.gz  https://i10git.cs.fau.de/api/v4/projects/walberla%2Fwalberla/repository/archive?sha=59c9b8b1
wget -LO cache/sources/e/ESPResSo/heffte-2.4.1.tar.gz       https://github.com/icl-utk-edu/heffte/archive/v2.4.1.tar.gz
wget -LO cache/sources/e/ESPResSo/kokkos-18b830e.tar.gz     https://github.com/kokkos/kokkos/archive/18b830e.tar.gz
wget -LO cache/sources/e/ESPResSo/Cabana-ebfaa51.tar.gz     https://github.com/ECP-copa/Cabana/archive/ebfaa51.tar.gz
wget -LO cache/sources/e/ESPResSo/highfive-0103467.tar.gz   https://github.com/highfive-devs/highfive/archive/0103467.tar.gz

time eb --force --include-easyblocks ./espresso.py ./ESPResSo-5.0-foss-2023b.eb --robot --prefix="${BUILDDIR}/cache" --skip-test-step
module use "${BUILDDIR}/cache/modules/all"
module load ESPResSo/4b9ccc4-foss-2023b
pypresso -c "import espressomd.version;print(espressomd.version.friendly())"
```

