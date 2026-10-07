"""Neutral simulations, step 3: reproducible SNP subset plans.

Sample k distinct neutral SNPs per technical replicate, where k is the tested gene's
variant count. Fish, locations and diversity metrics are handled by the Rmd.
Uses only Python's standard library; no bcftools commands are executed here.

The helper commands support step 2's Bash bcftools operations: they prepare
coordinate lists, check VCF records and record the selected subsets for the Rmd.
Run `plan --help` for sampling, or `--help` for those helper commands.
"""
import argparse
import csv
import gzip
import hashlib
import json
import random
import re
from pathlib import Path


def open_text(path):
    if str(path).lower().endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8-sig")


def fingerprint(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5_fingerprint(path):
    return fingerprint(path, "md5")


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def canonical_chrom(chrom):
    return re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)


def vcf_records(path):
    with open_text(path) as stream:
        for line in stream:
            if not line.startswith("#") and line.strip():
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 8:
                    raise ValueError(f"Malformed VCF record in {path}")
                yield fields, line


def vcf_sites(path):
    """Stream SNP identities without retaining genotype columns."""
    with open_text(path) as stream:
        for line in stream:
            if not line.startswith("#") and line.strip():
                fields = line.rstrip("\r\n").split("\t", 5)
                if len(fields) < 6:
                    raise ValueError(f"Malformed VCF record in {path}")
                yield fields[0], int(fields[1]), fields[3], fields[4]


def vcf_sample_header(path):
    with open_text(path) as stream:
        for line in stream:
            if line.startswith("#CHROM\t"):
                return line.rstrip("\r\n")
    raise ValueError(f"VCF sample header missing from {path}")


def is_biallelic_snp(ref, alt):
    return len(ref) == len(alt) == 1 and ref in "ACGT" and alt in "ACGT" and ref != alt


def guard_paths(args):
    inputs = {Path(path).resolve() for path in args.input}
    outputs = [Path(path).resolve() for path in args.output]
    if len(outputs) != len(set(outputs)) or inputs.intersection(outputs):
        raise ValueError("Neutral inputs and outputs must be distinct files; an output would overwrite an input or another output")


def read_markers(path, input_format):
    """Read converted tables or ready neutral coordinates, without table conversion."""
    markers = {}
    with open_text(path) as stream:
        rows = csv.DictReader(stream, delimiter="\t") if input_format == "normalized" else csv.reader(stream, delimiter="\t")
        if input_format == "normalized" and not {"CHROM", "POS"}.issubset(rows.fieldnames or []):
            raise ValueError("Converted neutral table requires CHROM and POS headers")
        for number, row in enumerate(rows, 2 if input_format == "normalized" else 1):
            if input_format == "normalized":
                chrom, position = (str(row.get(key) or "").strip() for key in ("CHROM", "POS"))
            else:
                if not row or not any(value.strip() for value in row):
                    continue
                if row[0].lstrip().startswith("#") or row[0].strip() in ("track", "browser") or row[0].startswith(("track ", "browser ")):
                    continue
                if input_format == "positions":
                    if len(row) != 2:
                        raise ValueError(f"Row {number}: positions requires exactly two tab-separated columns CHROM POS")
                    chrom, position = (value.strip() for value in row)
                else:
                    if len(row) < 3:
                        raise ValueError(f"Row {number}: BED requires CHROM START END")
                    chrom, start, position = (value.strip() for value in row[:3])
                    if not start.isdigit() or not position.isdigit() or int(position) != int(start) + 1:
                        raise ValueError(f"Row {number}: BED neutral sites must be single-base SNP intervals")
            if not chrom or not position.isdigit() or int(position) < 1:
                raise ValueError(f"Invalid CHROM/POS on row {number}")
            key = (chrom, int(position))
            markers[key] = None
    if not markers:
        raise ValueError("No neutral markers supplied")
    return markers


def pool_target_markers(args):
    """Resolve neutral coordinates to the cohort VCF's contig names."""
    header = Path(args.vcf_header).read_text(encoding="utf-8")
    contigs = re.findall(r"^##contig=<ID=([^,>]+)", header, re.MULTILINE)
    if not contigs and args.contigs:
        contigs = list(dict.fromkeys(Path(args.contigs).read_text(encoding="utf-8").splitlines()))
    aliases = {}
    for contig in contigs:
        key = canonical_chrom(contig)
        if key in aliases and aliases[key] != contig:
            raise ValueError(f"VCF contains ambiguous contigs {aliases[key]} and {contig}")
        aliases[key] = contig
    wanted = {}
    for chrom, pos in read_markers(args.markers, args.markers_format):
        key = (canonical_chrom(chrom), pos)
        if key in wanted:
            raise ValueError(f"Ambiguous marker position after chromosome normalization: {key}")
        wanted[key] = aliases.get(key[0])
    Path(args.target_markers_output).parent.mkdir(parents=True, exist_ok=True)
    mapped = 0
    with open(args.target_markers_output, "w", encoding="utf-8", newline="\n") as stream:
        for (_, pos), vcf_chrom in wanted.items():
            if vcf_chrom:
                stream.write(f"{vcf_chrom}\t{pos}\n")
                mapped += 1
    if not mapped:
        raise ValueError("No marker chromosome names match the VCF; check the marker coordinates and VCF contig names")
    print(f"Mapped {mapped}/{len(wanted)} neutral marker coordinates to VCF contigs")


def check_pool_records(args):
    """Keep biallelic SNPs and check that each selected position occurs once."""
    wanted = read_markers(args.target_markers, "positions")
    matched = set()
    with open_text(args.candidates) as source, open(args.validated_vcf, "w", encoding="utf-8", newline="\n") as dest:
        for line in source:
            if line.startswith("#"):
                dest.write(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise ValueError("Malformed VCF record")
            key = (fields[0], int(fields[1]))
            if key not in wanted:
                raise ValueError("bcftools returned a position outside the neutral target list")
            if not is_biallelic_snp(fields[3], fields[4]):
                continue
            if key in matched:
                raise ValueError(f"Multiple VCF records match neutral position {key}; supply a VCF with one SNP record per position")
            matched.add(key)
            dest.write(line)
    if not matched:
        raise ValueError("No neutral biallelic SNPs matched; check marker coordinates and the VCF records")
    print(f"Retained {len(matched)}/{len(wanted)} targeted positions as neutral biallelic SNPs")


def check_pool(args):
    """Validate a supplied neutral VCF without rewriting it."""
    source = Path(args.vcf)
    if not str(source).lower().endswith((".vcf.gz", ".vcf.bgz")):
        raise ValueError("Supplied neutral VCF must be a BGZF-compressed .vcf.gz or .vcf.bgz file")
    header = vcf_sample_header(source).split("\t")
    if len(header) < 10 or header[8] != "FORMAT" or any(not sample for sample in header[9:]):
        raise ValueError("Supplied neutral VCF requires a FORMAT column and named genotype samples")
    if len(set(header[9:])) != len(header[9:]):
        raise ValueError("Supplied neutral VCF has duplicate sample names")
    count, previous_chrom, previous_pos = 0, None, None
    chromosomes, aliases = set(), {}
    for count, (chrom, pos, ref, alt) in enumerate(vcf_sites(source), 1):
        if not chrom or pos < 1 or not is_biallelic_snp(ref, alt):
            raise ValueError(f"Supplied neutral VCF record {count} must be a biallelic SNP at a positive position")
        canonical = canonical_chrom(chrom)
        if canonical in aliases and aliases[canonical] != chrom:
            raise ValueError(f"Supplied neutral VCF contains ambiguous contigs {aliases[canonical]} and {chrom}")
        aliases[canonical] = chrom
        if chrom == previous_chrom:
            if pos == previous_pos:
                raise ValueError(f"Supplied neutral VCF has duplicate position {chrom}:{pos}")
            if pos < previous_pos:
                raise ValueError("Supplied neutral VCF must be sorted by position for indexed extraction")
        elif chrom in chromosomes:
            raise ValueError("Supplied neutral VCF must contain contiguous chromosome blocks for indexed extraction")
        chromosomes.add(chrom)
        previous_chrom, previous_pos = chrom, pos
    if not count:
        raise ValueError("Supplied neutral VCF contains no SNPs")
    print(f"Validated {count} supplied neutral biallelic SNPs for {len(header)-9} samples")


def plan(args):
    """Draw seeded, SNP-count-matched neutral subsets for one gene."""
    guard_paths(argparse.Namespace(input=[args.pool, args.gene_vcf], output=[args.output, args.metadata]))
    if args.technical_replicates < 1 or args.seed < 0:
        raise ValueError("technical_replicates must be positive and seed nonnegative")
    targets = list(vcf_sites(args.gene_vcf))
    if len({canonical_chrom(chrom) for chrom, _, _, _ in targets}) > 1:
        raise ValueError(f"Target VCF for {args.gene!r} contains multiple chromosomes; "
                         "provide the filtered VCF for the studied chromosome")
    if any(not is_biallelic_snp(ref, alt) for _, _, ref, alt in targets):
        raise ValueError("Neutral comparison requires biallelic SNPs in the target VCF")
    if len(set(targets)) != len(targets):
        raise ValueError("Gene VCF contains duplicate variant records")
    target_positions = {(canonical_chrom(t[0]), t[1]) for t in targets}
    available, pool_count = [], 0
    for pool_count, (chrom, pos, _, _) in enumerate(vcf_sites(args.pool), 1):
        if (canonical_chrom(chrom), pos) not in target_positions:
            available.append(pool_count)
    k = len(targets)
    status = "ok" if 0 < k <= len(available) else "no_target_variants" if not k else "insufficient_neutral_markers"
    gene_seed = int.from_bytes(hashlib.sha256(f"{args.seed}:{args.gene}".encode()).digest()[:8], "big")
    rng = random.Random(gene_seed)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["replicate", "marker_indices"])
        if status == "ok":
            for i in range(1, args.technical_replicates + 1):
                writer.writerow([i, ",".join(map(str, sorted(rng.sample(available, k))))])
    write_json(args.metadata, dict(gene_id=args.gene, status=status, n_variants=k, n_available=len(available),
        n_excluded_overlap=pool_count-len(available), pool_n_matched=pool_count, technical_replicates_requested=args.technical_replicates,
        seed=args.seed, gene_seed=str(gene_seed), sampling="without_replacement_within_replicate",
        pool_sha256=fingerprint(args.pool), gene_vcf_sha256=fingerprint(args.gene_vcf),
        pool_md5=md5_fingerprint(args.pool), gene_vcf_md5=md5_fingerprint(args.gene_vcf),
        plan_md5=md5_fingerprint(args.output), plan_sha256=fingerprint(args.output)))
    print(f"{args.gene}: {status}; target={k}, available={len(available)}")


def prepare_subsets(args):
    """Write the exact coordinate lists that step 2 passes to bcftools."""
    metadata = json.loads(Path(args.metadata).read_text(encoding="utf-8"))
    hashes = {"plan_md5": md5_fingerprint(args.plan), "pool_md5": md5_fingerprint(args.pool),
              "gene_vcf_md5": md5_fingerprint(args.gene_vcf)}
    for key, value in hashes.items():
        if metadata.get(key) != value:
            raise ValueError(f"{key} differs from plan metadata; regenerate neutral subsets before materializing")
    output_dir, work_dir = Path(args.output_dir), Path(args.work_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    protected = {Path(path).resolve() for path in (args.pool, args.gene_vcf, args.plan, args.metadata)}
    if args.pool_index:
        protected.add(Path(args.pool_index).resolve())
    if (output_dir / "manifest.json").resolve() in protected:
        raise ValueError("Prepared manifest would overwrite an input; choose a separate output directory")
    manifest = dict(status=metadata["status"], **hashes, n_variants=metadata["n_variants"],
                    n_technical_replicates=0, vcfs=[], replicate_files=[])
    groups = []
    if metadata["status"] == "ok":
        k = int(metadata["n_variants"])
        subsets = []
        with open(args.plan, encoding="utf-8-sig") as stream:
            for expected_rep, row in enumerate(csv.DictReader(stream, delimiter="\t"), 1):
                indices = [int(value) for value in row["marker_indices"].split(",")]
                if int(row["replicate"]) != expected_rep or len(indices) != k or len(set(indices)) != k or min(indices) < 1:
                    raise ValueError("Neutral subset plan must contain ordered replicates of exactly k distinct positive indices")
                subsets.append(indices)
        if len(subsets) != int(metadata["technical_replicates_requested"]):
            raise ValueError("Neutral subset plan does not contain the requested number of replicates")
        selected = set().union(*map(set, subsets))
        targets = {(canonical_chrom(chrom), pos) for chrom, pos, _, _ in vcf_sites(args.gene_vcf)}
        selected_sites, positions, count = {}, set(), 0
        for count, site in enumerate(vcf_sites(args.pool), 1):
            if count not in selected:
                continue
            chrom, pos, ref, alt = site
            if not is_biallelic_snp(ref, alt):
                raise ValueError(f"Selected neutral pool record {count} is not a biallelic SNP")
            if (canonical_chrom(chrom), pos) in targets:
                raise ValueError(f"Selected neutral pool record {count} overlaps the target gene")
            if (chrom, pos) in positions:
                raise ValueError(f"Selected neutral pool has multiple records at {chrom}:{pos}")
            positions.add((chrom, pos))
            selected_sites[count] = site
        if count != int(metadata["pool_n_matched"]) or selected != selected_sites.keys():
            raise ValueError("Neutral subset plan indices or pool record count disagree with the source pool")
        if args.pool_index:
            if not Path(args.pool_index).is_file():
                raise ValueError("Configured neutral pool index is missing; regenerate its cache CSI index")
        elif not any(Path(str(args.pool) + suffix).is_file() for suffix in (".csi", ".tbi")):
            raise ValueError("Indexed neutral pool required for materialization; regenerate its CSI index")
        groups = [
            (f"replicate{rep:04d}.vcf.gz", indices) for rep, indices in enumerate(subsets, 1)]
        for filename, _ in groups:
            if any((output_dir / (filename + suffix)).resolve() in protected for suffix in ("", ".csi")):
                raise ValueError("Prepared VCF output would overwrite a neutral source or plan input; choose a separate output directory")
        for filename, indices in groups:
            with open(work_dir / (filename + ".positions.tsv"), "w", encoding="utf-8", newline="\n") as stream:
                for index in sorted(indices):
                    chrom, pos, _, _ = selected_sites[index]
                    stream.write(f"{chrom}\t{pos}\n")
            manifest["vcfs"].append(dict(filename=filename, pool_indices=list(indices),
                                         expected_sites=[selected_sites[index] for index in indices]))
        manifest.update(n_technical_replicates=len(subsets), sample_header=vcf_sample_header(args.pool),
            replicate_files=[filename for filename, _ in groups])
    # Invalidate completion only after all input/output collision checks pass.
    (output_dir / "manifest.json").unlink(missing_ok=True)
    # Bash reads this protocol one filename per line, including on Windows.
    (work_dir / "groups.tsv").write_text(
        "".join(filename + "\n" for filename, _ in groups), encoding="utf-8", newline="\n")
    write_json(work_dir / "manifest.json", manifest)


def check_subset(args):
    """Match extracted SNP identities and sample order to a subset plan."""
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    entry = next(entry for entry in manifest["vcfs"] if entry["filename"] == args.filename)
    identities = {tuple(site): index for site, index in zip(entry["expected_sites"], entry["pool_indices"])}
    emitted = []
    for site in vcf_sites(args.vcf):
        if site not in identities:
            raise ValueError(f"Prepared VCF {args.filename} contains an unrequested variant {site}")
        emitted.append(identities[site])
    if len(emitted) != len(identities) or set(emitted) != set(identities.values()):
        raise ValueError(f"Prepared VCF {args.filename} has missing or duplicated selected SNPs")
    if vcf_sample_header(args.vcf) != manifest["sample_header"]:
        raise ValueError(f"Prepared VCF {args.filename} altered the sample order")
    entry["pool_indices"] = emitted
    write_json(args.manifest, manifest)


def finish_subsets(args):
    """Write Rmd-compatible provenance for completed subset VCFs."""
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    directory = Path(args.manifest).parent
    for entry in manifest["vcfs"]:
        path = directory / entry["filename"]
        if not Path(str(path) + ".csi").is_file():
            raise ValueError(f"Prepared VCF {path.name} has no CSI index")
        entry["md5"] = md5_fingerprint(path)
        entry.pop("expected_sites")
    manifest.pop("sample_header", None)
    if manifest["status"] == "ok":
        manifest.update(bcftools_version=args.bcftools_version, extraction="bcftools_view_indexed_regions_exact_POS")
    write_json(args.output, manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    def command(name, function, required):
        sub = commands.add_parser(name, help=function.__doc__)
        for flag in required.split():
            sub.add_argument("--" + flag, required=True)
        sub.set_defaults(func=function)
        return sub

    sub = command("plan", plan, "pool gene-vcf gene output metadata")
    sub.add_argument("--technical-replicates", type=int, default=100)
    sub.add_argument("--seed", type=int, default=42)
    sub = commands.add_parser("guard-paths", help="Check that outputs cannot overwrite inputs")
    sub.add_argument("--input", action="append", required=True)
    sub.add_argument("--output", action="append", required=True)
    sub.set_defaults(func=guard_paths)
    sub = command("pool-target-markers", pool_target_markers, "vcf-header markers target-markers-output")
    sub.add_argument("--markers-format", choices=("normalized", "positions", "bed"), default="normalized")
    sub.add_argument("--contigs", default="")
    command("check-pool-records", check_pool_records, "candidates target-markers validated-vcf")
    command("check-pool", check_pool, "vcf")
    sub = command("prepare-subsets", prepare_subsets, "pool gene-vcf plan metadata output-dir work-dir")
    sub.add_argument("--pool-index", default="")
    command("check-subset", check_subset, "vcf manifest filename")
    command("finish-subsets", finish_subsets, "manifest output bcftools-version")
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
