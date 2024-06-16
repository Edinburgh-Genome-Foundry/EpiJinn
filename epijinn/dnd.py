# Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of EpiJinn.
#
# EpiJinn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# EpiJinn is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

from .Methyl import Methylase


class Dnd(Methylase):
    pass


Dnd_EcoB7A = Dnd("Dnd_EcoB7A", "GAAC", 0, 3)
Dnd_Sli1326 = Dnd("Dnd_Sli1326", "GGCC", 0, 3)
Dnd_VciFF75 = Dnd("Dnd_VciFF75", "CCA", 0, 2)

DND = [Dnd_EcoB7A, Dnd_Sli1326, Dnd_VciFF75]
