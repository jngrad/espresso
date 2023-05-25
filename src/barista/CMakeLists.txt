#
# Copyright (C) 2023 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

add_library(espresso_barista SHARED)
add_library(espresso::barista ALIAS espresso_barista)
target_include_directories(espresso_barista PUBLIC include)

install(TARGETS espresso_barista
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)

add_subdirectory(src)

if(ESPRESSO_BUILD_TESTS)
  add_subdirectory(tests)
endif()
