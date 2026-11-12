# Copyright (C) 2026 Daniel Gamermann
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

from platformdirs import user_data_dir
from pathlib import Path

APP_NAME = "KEGGDB"
APP_AUTHOR = "DanielGamermann"

data_dir = Path(user_data_dir(APP_NAME, APP_AUTHOR))
data_dir.mkdir(parents=True, exist_ok=True)


keggPath = data_dir

obgs = ["compounds", "glycans", "reactions", "enzymes", 'modules', 'pathways', 'other']

for obg in obgs:
    path = keggPath / obg
    path.mkdir(parents=True, exist_ok=True)

