#
# Copyright 2025 Jean-Noël Grad
#
# This code is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
EasyBuild support for ESPResSo, implemented as an easyblock.

@author: Jean-Noël Grad (University of Stuttgart)
"""

import os
import re
from easybuild.easyblocks.generic.cmakemake import CMakeMake
from easybuild.tools.systemtools import get_cpu_architecture, get_cpu_features
from easybuild.tools.systemtools import X86_64


class EB_ESPResSo(CMakeMake):
    """Support for building and installing ESPResSo."""

    def _get_extracted_tarball_paths(self):
        extracted_paths = {}
        for src in self.src:
            name = src['name'].split('-', 1)[0]
            # process main software
            if name == 'espresso':
                extracted_paths['espresso'] = src['finalpath']
                continue
            # process dependencies
            tarball = src['name']
            if not tarball.endswith('.tar.gz'):
                raise ValueError(tarball + ' is not a tar.gz file')
            prefix = tarball.rsplit('.', 2)[0]
            matches = [x for x in os.listdir(src['finalpath']) if x.startswith(prefix)]
            if len(matches) == 0:
                raise RuntimeError(tarball + ' was not extracted')
            if len(matches) > 1:
                raise RuntimeError(tarball + ' matches multiple folders: ' + str(matches))
            extracted_paths[name] = os.path.join(src['finalpath'], matches[0])
        return extracted_paths

    def _patch_fetchcontent(self, extracted_paths):
        cmakelists_path = os.path.join(extracted_paths['espresso'], 'CMakeLists.txt')
        with open(cmakelists_path, 'r') as f:
            content = f.read()
        for name, local_uri in extracted_paths.items():
            if name == 'espresso':
                continue
            pattern = r'FetchContent_Declare\(\s*' + name + '\s+GIT_REPOSITORY\s+\S+\s+GIT_TAG\s+\S+(?=\s|\))'
            m = re.search(pattern, content)
            if m is None:
                raise RuntimeError(name + ' is not part of the ESPResSo FetchContent workflow')
            content = re.sub(pattern, 'FetchContent_Declare(' + name + ' URL ' + local_uri, content)
        with open(cmakelists_path, 'w') as f:
            f.write(content)

    def configure_step(self):
        # patch FetchContent to avoid downloading dependencies
        extracted_paths = self._get_extracted_tarball_paths()
        self._patch_fetchcontent(extracted_paths)

        configopts = self.cfg.get('configopts', '')
        dependencies = self.cfg.get('dependencies', [])

        cpu_features = get_cpu_features()
        with_cuda = any(x.get('name', '') == 'CUDA' for x in dependencies)
        with_pfft = any(x.get('name', '') == 'PFFT' for x in dependencies)
        with_hdf5 = any(x.get('name', '') == 'HDF5' for x in dependencies)
        with_gsl = any(x.get('name', '') == 'GSL' for x in dependencies)

        if with_cuda:
            configopts += ' -DESPRESSO_BUILD_WITH_CUDA=ON'
        else:
            configopts += ' -DESPRESSO_BUILD_WITH_CUDA=OFF'
        if with_hdf5:
            configopts += ' -DESPRESSO_BUILD_WITH_HDF5=ON'
        else:
            configopts += ' -DESPRESSO_BUILD_WITH_HDF5=OFF'
        if with_gsl:
            configopts += ' -DESPRESSO_BUILD_WITH_GSL=ON'
        else:
            configopts += ' -DESPRESSO_BUILD_WITH_GSL=OFF'

        configopts += ' -DESPRESSO_BUILD_WITH_WALBERLA=ON'
        if with_pfft:
            configopts += ' -DESPRESSO_BUILD_WITH_WALBERLA_FFT=ON'
        else:
            configopts += ' -DESPRESSO_BUILD_WITH_WALBERLA_FFT=OFF'
        if get_cpu_architecture() == X86_64 and 'avx2' in cpu_features:
            configopts += ' -DESPRESSO_BUILD_WITH_WALBERLA_AVX=ON'

        self.cfg['configopts'] = configopts

        return super(EB_ESPResSo, self).configure_step()
