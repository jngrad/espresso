#!/usr/bin/env sh
#
# Copyright (C) 2017-2025 The ESPResSo project
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

[ "${#}" -eq 1 ] || exit 1

GIT_COMMIT=$(git rev-parse HEAD)
TARGET_URL="https://gitlab.icp.uni-stuttgart.de/jgrad/espresso/-/pipelines/${CI_PIPELINE_ID}"
STATUS="${1}"
if [ "${STATUS}" = "pending" ]; then DESCRIPTION="The CI pipeline has started"; fi
if [ "${STATUS}" = "failure" ]; then DESCRIPTION="Failing"; fi
if [ "${STATUS}" = "success" ]; then DESCRIPTION="Successful"; fi
curl -L "https://api.github.com/repos/jngrad/espresso/statuses/${GIT_COMMIT}" \
     -H "Authorization: token ${GITHUB_TOKEN}" \
     -H "Accept: application/vnd.github+json" \
     -H "X-GitHub-Api-Version: 2022-11-28" \
     -X POST \
     -d '{"state":"'${STATUS}'", "context":"ICP GitLab CI", "description":"'${DESCRIPTION}'", "target_url":"'${TARGET_URL}'"}'
