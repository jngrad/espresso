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

from .script_interface import script_interface_register, ScriptInterfaceHelper


@script_interface_register
class FancyBrew(ScriptInterfaceHelper):
    """
    Brew a fancy beverage.

    Parameters
    ----------
    base : :obj:`str`
        Can be one of 'espresso', 'ristretto', 'lungo'.
    toppings : :obj:`list`
        Can be any combination of 'milk', 'sugar', 'cream'.

    Methods
    -------
    brew_coffee()
        Brew and serve the beverage with the currently active flavoring options.

        Returns
        -------
        : :obj:`str`
            A description of the beverage.

    set_toppings()
        Set the flavoring options.

        Parameters
        ----------
        toppings : :obj:`list`
            The flavoring options.

    run_cleaning_cycle()
        Clean the brewing chamber and if necessary the milk steamer.
        In a MPI-parallel simulation, the volume of water grows linearly
        with the number of cores.

        Returns
        -------
        : :obj:`float`
            The volume of water used in liters.

    """
    _so_name = "Barista::FancyBrew"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "brew_coffee",
        "set_toppings",
        "run_cleaning_cycle",
    )
