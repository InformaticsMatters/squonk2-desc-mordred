# Copyright 2024 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Mordred featurization that reproduces jaqpotpy 1.0.x, for the legacy
environment (original mordred 1.2.0, rdkit 2023.9.6, numpy 1.26).

jaqpotpy built ``Calculator(descriptors, ignore_3D=True)`` and, per molecule,
kept each descriptor value that was a plain number and turned everything else
(mordred's Missing/Error objects) into NaN. It also featurised the
SMILES-canonical molecule (strip explicit Hs, canonicalise, re-parse). A NaN
cell is written empty, which the predictor fills with the model's training fill
(0 for the legacy models) — the same end result jaqpotpy produced.

This matches the predictor's proven Mordred featurizer. mordredcommunity (what
`main` uses) changes some values, which is why reproduction needs this branch.
"""

import numpy
from mordred import Calculator, descriptors
from rdkit import Chem

# the predictor's featurizer name, written into the manifest's featurizer.name
MORDRED = "mordred"

# ignore_3D is the one value-affecting knob jaqpotpy set (and the predictor
# checks); the legacy models were all 2D.
_CALCULATOR = Calculator(descriptors, ignore_3D=True)
_NAMES = [str(d) for d in _CALCULATOR.descriptors]

_NUMERIC = (int, float, numpy.integer, numpy.floating, numpy.bool_)


def descriptor_names():
    """The 1613 mordred descriptor column names, in order."""
    return list(_NAMES)


def mordred_params():
    """The value-affecting mordred knobs, as the manifest/predictor want them."""
    return {"ignore_3D": True}


def normalize(mol):
    """Bring a molecule to what parsing its canonical SMILES would give, as
    jaqpotpy did before featurising. Raises if it can't be processed, which the
    caller turns into a flagged (empty) output row."""
    normalized = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.RemoveHs(mol)))
    if normalized is None or normalized.GetNumAtoms() == 0:
        raise ValueError("molecule could not be normalised")
    return normalized


def mordred_values(mol):
    """Descriptor values for one molecule, in ``descriptor_names`` order.

    Each cell is the numeric value if mordred returned a finite number, else
    None (written as an empty cell and filled by the predictor) — this is how
    jaqpotpy treated mordred's Missing/Error results and non-finite values.
    """
    result = _CALCULATOR(normalize(mol))
    out = []
    for value in result.values():
        if isinstance(value, _NUMERIC):
            number = float(value)
            out.append(number if numpy.isfinite(number) else None)
        else:
            out.append(None)
    return out
