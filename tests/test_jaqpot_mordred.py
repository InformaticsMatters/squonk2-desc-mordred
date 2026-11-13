"""Reference-value tests for the legacy mordred featurizer.

These freeze the output of ``jaqpot_mordred`` under the branch's pinned stack
(original mordred 1.2.0, rdkit 2023.9.6, numpy 1.26), so the build fails loudly
if anything drifts — in particular a switch to mordredcommunity, which changes
some values. The golden values were generated once with that stack, which
matches the predictor's proven Mordred featurizer.
"""

import pytest
from rdkit import Chem

import jaqpot_mordred as jm


def _mol(smiles):
    return Chem.MolFromSmiles(smiles)


def test_descriptor_count_and_order():
    names = jm.descriptor_names()
    assert len(names) == 1613
    assert names[:3] == ["ABC", "ABCGG", "nAcid"]


def test_params_match_the_manifest():
    assert jm.mordred_params() == {"ignore_3D": True}


def test_ethanol_values():
    values = dict(zip(jm.descriptor_names(), jm.mordred_values(_mol("CCO"))))
    assert values["MW"] == pytest.approx(46.041864812, abs=1e-6)
    assert values["nAtom"] == 9.0
    # ABC needs a ring; mordred returns a missing value for acyclic ethanol,
    # which we write as None (empty cell, filled by the predictor)
    assert values["ABC"] is None


def test_aspirin_values():
    values = dict(
        zip(jm.descriptor_names(), jm.mordred_values(_mol("CC(=O)Oc1ccccc1C(=O)O")))
    )
    assert values["nAromAtom"] == 6.0
    assert values["nC"] == 9.0
    assert values["nBase"] == 0.0
    # original mordred 1.2.0 returns a missing value for ABC/ABCGG here; this
    # is exactly the kind of difference that makes the legacy branch necessary
    assert values["ABC"] is None


def test_explicit_hs_match_plain_smiles():
    assert jm.mordred_values(Chem.AddHs(_mol("CCO"))) == jm.mordred_values(_mol("CCO"))


def test_normalize_raises_on_unusable_molecule():
    with pytest.raises(ValueError):
        jm.normalize(Chem.Mol())
