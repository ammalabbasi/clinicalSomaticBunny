#!/bin/bash
# Install PCGR 2.3.2 (conda, the officially supported non-root route) and its public GRCh38/GRCh37 data.
# Public data only: the PCGR reference bundle (insilico.hpc.uio.no) and the Ensembl VEP 115 indexed cache.
#
#   bash Modules/PCGR/install_pcgr.sh <dest> envs              # conda envs <dest>/conda/pcgr + <dest>/conda/pcgrr
#   bash Modules/PCGR/install_pcgr.sh <dest> refdata grch38    # <dest>/refdata/20260620/data/grch38   (~7.7 GB download)
#   bash Modules/PCGR/install_pcgr.sh <dest> vep GRCh38        # <dest>/vep_cache/homo_sapiens/115_GRCh38 (~25.3 GB download, tarball kept until you delete it)
#
# Then run SomaticBunny with
#   --pcgr_env <dest>/conda/pcgr --pcgr_refdata <dest>/refdata/20260620 --pcgr_vep_dir <dest>/vep_cache
#
# The two envs MUST stay siblings named pcgr and pcgrr: pcgr finds the R env as $(dirname $CONDA_PREFIX)/pcgrr.
# The bundle version (20260620) and VEP release (115) are pinned by PCGR 2.3.2; do not mix versions.
set -euo pipefail
# A user-site pandas (~/.local/lib/python3.12) shadows the env and breaks pcgr ("No module named pytz").
export PYTHONNOUSERSITE=1
unset PYTHONPATH

DEST=${1:?usage: install_pcgr.sh <dest> envs|refdata|vep [assembly]}
WHAT=${2:?usage: install_pcgr.sh <dest> envs|refdata|vep [assembly]}
PCGR_VERSION=2.3.2
BUNDLE_VERSION=20260620
VEP_VERSION=115
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

case "$WHAT" in
  envs)
    mkdir -p "$DEST/conda"
    # Lock files copied from github.com/sigven/pcgr/tree/v${PCGR_VERSION}/conda/env/lock (explicit, fully pinned)
    conda create -y --prefix "$DEST/conda/pcgr"  --file "$HERE/yml/pcgr/pcgr-linux-64.lock"
    conda create -y --prefix "$DEST/conda/pcgrr" --file "$HERE/yml/pcgr/pcgrr-linux-64.lock"
    "$DEST/conda/pcgr/bin/pcgr" --version
    ;;
  refdata)
    GENOME=${3:-grch38}
    mkdir -p "$DEST/refdata/$BUNDLE_VERSION"
    cd "$DEST/refdata/$BUNDLE_VERSION"
    # Streamed: the .tgz is never stored, tar fails on a truncated stream.
    curl -fL --retry 5 "https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.${BUNDLE_VERSION}.${GENOME}.tgz" | tar xzf -
    cat "data/${GENOME}/.PCGR_BUNDLE_VERSION"
    ;;
  vep)
    GENOME=${3:-GRCh38}
    mkdir -p "$DEST/vep_cache"
    cd "$DEST/vep_cache"
    # Download to a file first (resumable with -C -; a single stream from ftp.ensembl.org can drop to <1 MB/s,
    # so a parallel range download is fine too), check Ensembl's BSD sum, then extract.
    BASE="https://ftp.ensembl.org/pub/release-${VEP_VERSION}/variation/indexed_vep_cache"
    TGZ="homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz"
    curl -fL --retry 5 -C - -o "$TGZ" "$BASE/$TGZ"
    curl -fsL "$BASE/CHECKSUMS" | grep "$TGZ"
    sum "$TGZ"    # must equal the two numbers printed above
    tar xzf "$TGZ"
    ls -d "homo_sapiens/${VEP_VERSION}_${GENOME}"
    echo "Extracted. Remove $DEST/vep_cache/$TGZ after checking the sums match."
    ;;
  *)
    echo "unknown step $WHAT (envs | refdata | vep)" >&2; exit 1 ;;
esac
