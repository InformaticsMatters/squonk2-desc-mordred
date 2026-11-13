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


# Generates Mordred descriptors for molecules.

import argparse
import json
import os
import sys
import time

import mordred
import numpy
import rdkit
import rdkit_utils
from dm_job_utilities.cli import ProgressReporter, add_reporting_args
from dm_job_utilities.dm_log import DmLog
from dm_job_utilities.utils import expand_path, is_type, log, read_delimiter
from rdkit import Chem

import jaqpot_mordred as jm

MANIFEST_SCHEMA_VERSION = 1


def manifest_path(outfile):
    """The sidecar manifest path for a descriptor file, matching how the
    predictor looks for it: the final extension (and a trailing ``.gz``) is
    replaced by ``.manifest.json`` (``foo.smi`` -> ``foo.manifest.json``,
    ``foo.smi.gz`` -> ``foo.manifest.json``)."""
    stem = outfile
    if stem.endswith(".gz"):
        stem = stem[: -len(".gz")]
    stem, _ = os.path.splitext(stem)
    return stem + ".manifest.json"


def write_manifest(outfile, mode, id_column, columns):
    """Write the sidecar manifest describing the mordred descriptor file, so
    the predictor can check the file fits the model."""
    is_sdf = outfile.endswith(".sdf") or outfile.endswith(".sdf.gz")

    id_is_int, _ = is_type(id_column, int)
    if id_column is None:
        id_name = "_Name" if is_sdf else rdkit_utils.SMILES_COL_NAME
    elif id_is_int or id_column == "_Name":
        id_name = rdkit_utils.ID_COL_NAME
    else:
        id_name = id_column
    smiles_column = None if is_sdf else rdkit_utils.SMILES_COL_NAME

    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "featurizer": {"name": jm.MORDRED, "params": jm.mordred_params()},
        "columns": list(columns),
        "id_column": id_name,
        "smiles_column": smiles_column,
        "fragment_method": mode,
        "library_versions": {
            "rdkit": rdkit.__version__,
            "numpy": numpy.__version__,
            "mordred": mordred.__version__,
        },
        # the job leaves uncomputable descriptors empty; the predictor fills
        # them with the model's training fill
        "missing_value": None,
    }
    path = manifest_path(outfile)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    DmLog.emit_event("Wrote descriptor manifest to {}".format(path))


def run(
    input_file,
    outfile,
    mode="hac",
    include_3d=False,
    delimiter=None,
    id_column=None,
    mol_column=None,
    omit_fields=False,
    read_header=False,
    write_header=False,
    read_records=50,
    interval=1000,
):

    if include_3d:
        if not (input_file.endswith(".sdf") or input_file.endswith(".sdf.gz")):
            DmLog.emit_event(
                "If calculating 3D descriptors input must be a SDF with 3D molecules"
            )
            exit(1)
        if not rdkit_utils.check_molecules_are_3d(input_file):
            DmLog.emit_event("Molecules do not seem to have 3D coordinates")
            exit(1)

    expand_path(outfile)

    count = 0
    errors = 0

    t0 = time.time()

    # setup the reader
    reader = rdkit_utils.create_reader(
        input_file,
        id_column=id_column,
        mol_column=mol_column,
        read_records=read_records,
        read_header=read_header,
        delimiter=delimiter,
    )
    extra_field_names = reader.get_extra_field_names()

    calc_field_names = jm.descriptor_names()
    DmLog.emit_event("Calculating {} descriptors".format(len(calc_field_names)))

    # setup the writer
    writer = rdkit_utils.create_writer(
        outfile,
        extra_field_names=extra_field_names,
        calc_prop_names=calc_field_names,
        delimiter=delimiter,
        id_column=id_column,
        mol_column=mol_column,
    )

    id_col_type, id_col_value = is_type(id_column, int)
    reporter = ProgressReporter(interval)

    # read the input records and write the output
    while True:
        t = reader.read()
        # break if no more data to read
        if not t:
            break
        mol, smi, mol_id, props = t

        if count == 0 and write_header:
            headers = rdkit_utils.generate_headers(
                id_col_type,
                id_col_value,
                reader.get_mol_field_name(),
                reader.field_names,
                calc_field_names,
                omit_fields,
            )

            writer.write_header(headers)

        count += 1

        reporter.report(count)
        if count % 50000 == 0:
            # Emit a 'total' cost, replacing all prior costs
            DmLog.emit_cost(count)

        if not mol:
            errors += 1
            DmLog.emit_event("Failed to read molecule for record", count)
            # keep the output lined up with the input: a flagged (empty) row
            # the predictor reports as a failed molecule
            _write_flagged_row(writer, smi, mol, mol_id, props, calc_field_names)
            continue

        if omit_fields:
            for name in mol.GetPropNames():
                mol.ClearProp(name)
            props = []

        try:
            if mode == "none":
                biggest = mol
                cann_smi = smi
            else:
                biggest = rdkit_utils.fragment(mol, mode)
                cann_smi = Chem.MolToSmiles(biggest)

            values = jm.mordred_values(biggest)

        except KeyboardInterrupt:
            log("Interrupted")
            sys.exit(0)

        except Exception:  # pylint: disable=broad-except
            errors += 1
            DmLog.emit_event("Failed to process record", count)
            _write_flagged_row(writer, smi, mol, mol_id, props, calc_field_names)
            continue

        writer.write(cann_smi, biggest, mol_id, props, values)

    writer.close()
    reader.close()

    write_manifest(outfile, mode, id_column, calc_field_names)

    t1 = time.time()
    log("Processing took {} secs".format(round(t1 - t0)))
    DmLog.emit_event("Processed {} records ({} failed)".format(count, errors))


def _write_flagged_row(writer, smi, mol, mol_id, props, columns):
    """Write a row with empty descriptor values, so a molecule that could not
    be read or featurised still appears in the output (lined up with the input)
    and is visibly a failure. Best-effort: if the row can't be represented it
    is dropped rather than crashing the run."""
    empty = [None] * len(columns)
    try:
        writer.write(smi, mol, mol_id, props, empty)
    except Exception as exc:  # pylint: disable=broad-except
        DmLog.emit_event("Could not write a flagged row: {}".format(exc))


def main():

    # Run using conda env created from environment-im-mordred.yaml
    #   or docker environment created from Dockerfile-mordred

    # Examples:
    #   python -m im_mordred.descriptor_generator -i data/10.smi -o descriptors.smi -d tab
    #   python -m im_mordred.descriptor_generator -i data/10+H.smi -o descriptors.smi -d tab --id_column 1 --read-header --writeheader

    # ----- command line args definitions ---------------------------------------------

    parser = argparse.ArgumentParser(description="Mordred 2D descriptors")
    # to pass tab as the delimiter specify it as $'\t' or use one of the
    # symbolic names 'comma', 'tab', 'space' or 'pipe'
    rdkit_utils.add_common_molecule_io_args(parser, output_default="descriptors2d.smi")
    add_reporting_args(parser)

    parser.add_argument(
        "-f",
        "--fragment-method",
        choices=["hac", "mw", "none"],
        default="hac",
        help="Strategy for picking largest fragment (mw or hac or none",
    )

    # parser.add_argument(
    #     "-3",
    #     "--include-3d",
    #     action="store_true",
    #     help="Include 3D descriptors (requires 3D molecules in SDF file)",
    # )

    args = parser.parse_args()
    DmLog.emit_event("descriptor_calc: ", args)

    delimiter = read_delimiter(args.delimiter)

    run(
        args.input,
        args.output,
        mode=args.fragment_method,
        # include_3d=args.include_3d,
        omit_fields=args.omit_fields,
        delimiter=delimiter,
        id_column=args.id_column,
        mol_column=args.mol_column,
        read_header=args.read_header,
        write_header=args.write_header,
        read_records=args.read_records,
        interval=args.interval,
    )


if __name__ == "__main__":
    main()
