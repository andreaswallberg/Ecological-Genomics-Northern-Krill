"""Prune VCF file by linkage disequilibrium

Prune a VCF file from variants that are in linkage. The default output is a
zarr directory that can be converted to vcf output for downstream processing
(e.g., as input to PLINK pca). Alternatively, bed output can be produced
for region-based filtering of vcf files.

This script is loosely based on the implementation in scikit-allel described
here: https://alimanfoo.github.io/2015/09/28/fast-pca.html

EXAMPLES

python ld_prune.py input.vcf.gz
python ld_prune.py --window-size

"""
import logging
import re
import tempfile
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
import sgkit as sg
import xarray as xr
from sgkit.io import vcf as sgvcf
from sgkit.typing import PathType

logger = logging.getLogger(__name__)


def convert_vcf_to_zarr(path: PathType, tmpdir: Optional[PathType] = None) -> PathType:
    """Convert vcf file to zarr"""
    output = Path(f"{path}.zarr")
    if tmpdir is not None:
        if isinstance(tmpdir, str):
            tmpdir = Path(tmpdir)
        elif isinstance(tmpdir, tempfile.TemporaryDirectory):
            tmpdir = Path(tmpdir.name)
        output = tmpdir.joinpath(output.name)
    logger.info(f"converting {path} to {output}...")
    sgvcf.vcf_to_zarr(path, output)
    return output


def to_bed(data: xr.Dataset, output: PathType) -> None:
    """Convert xarray dataset to bed file"""
    try:
        data_df = pd.DataFrame(
            [
                [data.contigs[i] for i in data.variant_contig.values],
                [str(i - 1) for i in data.variant_position.values],
                [str(i) for i in data.variant_position.values],
            ]
        ).T
    except Exception as exc:
        print(exc)
        raise
    print(f"converting dataset {repr(data)} to {output}")
    data_df.to_csv(output, index=False, sep="\t", header=False)


def run_ld_prune(vcf, **kwargs):
    """Run ld pruning on vcf data"""
    tmpdir = kwargs["tmpdir"]
    if tmpdir is None:
        tmpdir = tempfile.TemporaryDirectory()
    for vcfdata in vcf:
        zarrdata = convert_vcf_to_zarr(vcfdata, tmpdir=tmpdir)
        _ld_prune(zarrdata, **kwargs)


def _ld_prune(
    zarrdata,
    *,
    subsample=None,
    subsample_fraction=None,
    window_size=500,
    window_step=250,
    threshold=0.1,
    n_iter=5,
    output_suffix=None,
    as_bed=False,
    **kwargs,
):
    logger.info(f"loading dataset {zarrdata}...")
    data = sg.load_dataset(zarrdata)
    nvars = None
    if subsample is not None:
        nvars = subsample
    if subsample_fraction is not None:
        nvars = min(len(data.variants), int(len(data.variants) * subsample_fraction))
    if nvars is not None:
        vidx = sorted(np.random.choice(len(data.variants), nvars, replace=False))
        print(f"subsetting dataset to {nvars} variants...")
        data = data.isel(variants=sorted(vidx))

    # For rechunking
    original_chunk_size = data.chunks["variants"][0]
    data["dosage"] = data["call_genotype"].sum(dim="ploidy")
    try:
        data = sg.window_by_variant(data, size=window_size, step=window_step)
    except Exception as exc:
        print(exc)
        raise

    print("pruning by ld...")
    for i in range(n_iter):
        n_start = len(data.variants)
        data = sg.ld_prune(data, threshold=threshold)
        n_retain = len(data.variants)
        n_remove = n_start - n_retain
        print(
            f"ld_prune_sgkit: iteration {i + 1}; retaining {n_retain}, "
            f"removing {n_remove}"
        )
        if n_remove == 0:
            print(f"ld pruning has converged after {i+1} iterations; quitting")
            break
        data = data.drop_vars(["window_contig", "window_start", "window_stop"])
        data = sg.window_by_variant(data, size=window_size, step=window_step)

    # Save dataset
    re_zarr = re.compile("(.zarr)$")
    output = re_zarr.sub(f"{output_suffix}.zarr", str(zarrdata))
    data = data.chunk(original_chunk_size)

    try:
        data["variant_id"] = data.variant_id.astype(str)
        data["variant_allele"] = data.variant_allele.astype(str)
        data["sample_id"] = data.sample_id.astype(str)
        sg.save_dataset(data, output, mode="w")
    except Exception as exc:
        print(exc)
        raise
    if as_bed:
        bedout = re_zarr.sub(".bed", output)
        to_bed(data, bedout)


@click.command(help=__doc__, context_settings={"show_default": True})
@click.argument("vcf", nargs=-1)
@click.option(
    "--window-size",
    help=(
        "the window size, measured in number of snps, "
        "to use for identifying ld blocks"
    ),
    default=500,
    type=click.IntRange(
        1,
    ),
)
@click.option(
    "--window-step",
    help=(
        "the window step size, measured in number of snps, "
        "to use for identifying ld blocks"
    ),
    default=250,
    type=click.IntRange(
        1,
    ),
)
@click.option(
    "--threshold",
    help="maximum value of r^2 to include variants",
    default=0.1,
    type=click.FloatRange(min=0.0, max=1.0),
)
@click.option(
    "--n-iter",
    help="number of pruning iterations",
    default=5,
    type=click.IntRange(
        1,
    ),
)
@click.option(
    "--subsample",
    "-s",
    help="subsample raw input to this number of sites",
    type=click.IntRange(
        1,
    ),
)
@click.option(
    "--subsample-fraction",
    "-f",
    help="subsample raw input to this fraction of sites",
    type=click.FloatRange(min=0.0, max=1.0),
)
@click.option(
    "--exclude", "-e", help="list of samples to exclude", multiple=True, type=str
)
@click.option("--as-bed", help="save output as bed", is_flag=True)
@click.option("--output-file", "-o", help="Output file name", type=click.Path())
@click.option(
    "--output-suffix", help="Output file suffix", default=".ld_prune", type=str
)
@click.option(
    "--tmpdir", help="temporary directory to write to", default=None, type=click.Path()
)
def main(vcf, **kwargs):
    """Prune variants based on LD."""
    if len(vcf) == 0:
        logger.warning("No input vcf files supplied")
    run_ld_prune(vcf, **kwargs)


if __name__ == "__main__":
    main()
