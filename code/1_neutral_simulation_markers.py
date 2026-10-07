"""Step 1 (optional): convert a marker table into neutral simulation coordinates.

Read an annotated XLSX, TSV or CSV table, select neutral rows, and write a
CHROM/POS TSV. This script only converts marker inputs; it does not read VCFs,
call bcftools or sample Monte Carlo replicates.

The neutral workflow bypasses this conversion when given a ready two-column
CHROM/POS position list or a single-base BED file. Standalone conversion of
those formats is also supported here. Table and position-list coordinates are
1-based; BED coordinates are 0-based, half-open. Users supply coordinates in
the cohort VCF's assembly; no liftover or allele conversion is performed.

Uses only the Python standard library; no package installation is required.
XLSX cells are read from their saved values, without Excel or openpyxl.
"""

import argparse
import csv
import gzip
import posixpath
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).lower().endswith((".gz", ".bgz")) else open(path, encoding="utf-8-sig")


def xlsx_rows(path, sheet):
    """Read cached cell values from a selected XLSX worksheet (no Excel needed)."""
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as book:
        strings = []
        if "xl/sharedStrings.xml" in book.namelist():
            strings = ["".join(e.itertext()) for e in ET.fromstring(book.read("xl/sharedStrings.xml")).findall("s:si", ns)]
        sheets = ET.fromstring(book.read("xl/workbook.xml")).findall("s:sheets/s:sheet", ns)
        chosen = next((s for s in sheets if s.get("name") == sheet), None)
        if chosen is None:
            raise ValueError(f"Worksheet {sheet!r} missing; available: {[s.get('name') for s in sheets]}")
        rid = chosen.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
        rels = ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
        target = next(r.get("Target") for r in rels if r.get("Id") == rid)
        member = target.lstrip("/") if target.startswith("/") else posixpath.normpath("xl/" + target)
        root = ET.fromstring(book.read(member))
        rows = {}
        for row in root.findall("s:sheetData/s:row", ns):
            values = {}
            for cell in row.findall("s:c", ns):
                column = 0
                for c in re.match(r"[A-Z]+", cell.get("r")).group():
                    column = column * 26 + ord(c) - 64
                node = cell.find("s:v", ns)
                value = node.text if node is not None else ""
                if cell.get("t") == "s":
                    value = strings[int(value)]
                elif cell.get("t") == "inlineStr":
                    value = "".join(t.text or "" for t in cell.findall("s:is//s:t", ns))
                values[column - 1] = value
            rows[int(row.get("r"))] = [values.get(i, "") for i in range(max(values, default=-1) + 1)]
        return [rows.get(i, []) for i in range(1, max(rows, default=0) + 1)]


def column_index(selector, header):
    """Accept a 1-based column number or a case-insensitive header name."""
    selector = str(selector).strip()
    if selector.isdigit():
        index = int(selector) - 1
        if not 0 <= index < len(header):
            raise ValueError(f"Column {selector} is outside the {len(header)}-column header")
        return index
    matches = [i for i, name in enumerate(header) if str(name).strip().casefold() == selector.casefold()]
    if len(matches) != 1:
        raise ValueError(f"Expected one column {selector!r}; found headers {header}")
    return matches[0]


def resolve_header_row(rows, setting, chrom_column, pos_column):
    """Find one coordinate header, or use an explicit 1-based row number."""
    setting = str(setting).strip().casefold()
    if setting != "auto":
        if not setting.isdigit() or int(setting) < 1:
            raise ValueError("header-row must be 'auto' or a positive 1-based row number")
        number = int(setting)
        if number > len(rows):
            raise ValueError("header-row is beyond the end of the input")
        return number

    aliases = {"chrom": {"chr", "chrom", "chromosome"},
               "pos": {"pos", "position", "pos1", "position1"}}
    candidates = []
    for number, row in enumerate(rows, 1):
        for field, selector in (("chrom", chrom_column), ("pos", pos_column)):
            try:
                index = column_index(selector, row)
            except ValueError:
                break
            # Named selectors already require an exact case-insensitive match.
            # Numeric selectors require a complete header label, not a substring
            # in preamble text. POS 1 / POS_1 are accepted as POS1.
            if str(selector).strip().isdigit():
                label = re.sub(r"[\s_.()-]+", "", str(row[index]).strip().lstrip("#").casefold())
                if label not in aliases[field]:
                    break
        else:
            candidates.append(number)
    if len(candidates) != 1:
        detail = f"found candidate rows {candidates}" if candidates else "found no matching coordinate header"
        raise ValueError(f"Cannot auto-detect a unique header: {detail}. "
                         "Use Chr/CHROM/Chromosome and Pos/Position/Pos 1 in the configured columns, "
                         "or set header_row (--header-row) to an explicit row number.")
    return candidates[0]


def normalize(args):
    if Path(args.source).resolve() == Path(args.output).resolve():
        raise ValueError("Neutral marker output would overwrite its source; choose a separate output file")
    source_name = str(args.source).lower().removesuffix(".gz")
    input_format = args.input_format
    if input_format == "auto":
        input_format = "bed" if source_name.endswith(".bed") else "table"
    if source_name.endswith(".xlsx"):
        if input_format != "table":
            raise ValueError("XLSX input must use table format")
        rows = xlsx_rows(args.source, args.sheet)
    else:
        with open_text(args.source) as stream:
            rows = list(csv.reader(stream, delimiter="," if source_name.endswith(".csv") else "\t"))
    columns = None
    if input_format == "table":
        header_row = resolve_header_row(rows, args.header_row, args.chrom_column, args.pos_column)
        header = [str(x).strip() for x in rows[header_row - 1]]
        columns = {
            "chrom": column_index(args.chrom_column, header),
            "pos": column_index(args.pos_column, header),
        }
        if not args.neutral_only:
            columns["label"] = column_index(args.marker_label_column or "TYPE", header)
        print(f"Using table header at row {header_row}")
        numbered_rows = enumerate(rows[header_row:], header_row + 1)
    else:
        numbered_rows = enumerate(rows, 1)
    markers = {}
    for number, row in numbered_rows:
        if not any(str(x).strip() for x in row):
            continue
        if input_format != "table" and (row[0].lstrip().startswith("#") or row[0].strip() in ("track", "browser") or row[0].startswith(("track ", "browser "))):
            continue
        if input_format == "table":
            if len(row) <= max(columns.values()):
                raise ValueError(f"Row {number} is missing a configured column")
            if "label" in columns and str(row[columns["label"]]).strip().casefold() != args.neutral_marker_label.strip().casefold():
                continue
            chrom = str(row[columns["chrom"]]).strip()
            pos_raw = str(row[columns["pos"]]).strip()
        elif input_format == "positions":
            if len(row) != 2:
                raise ValueError(f"Row {number}: positions format requires exactly two tab-separated columns: CHROM POS")
            chrom, pos_raw = (str(x).strip() for x in row)
        else:
            if len(row) < 3:
                raise ValueError(f"Row {number}: BED requires CHROM START END; use --input-format positions for a two-column 1-based list")
            chrom, start, end = (str(x).strip() for x in row[:3])
            if not start.isdigit() or not end.isdigit() or int(end) != int(start) + 1:
                raise ValueError(f"Row {number}: neutral BED entries must be single-base SNP intervals [start, start+1)")
            pos_raw = end
        if not chrom or not re.fullmatch(r"[0-9]+", pos_raw) or int(pos_raw) < 1:
            raise ValueError(f"Invalid CHROM/POS on row {number}: {chrom!r}, {pos_raw!r}")
        key = (chrom, int(pos_raw))
        markers[key] = None
    if not markers:
        raise ValueError("No neutral markers selected")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as stream:
        out = csv.writer(stream, delimiter="\t", lineterminator="\n")
        out.writerow(["CHROM", "POS"])
        for chrom, pos in markers:
            out.writerow([chrom, pos])
    print(f"Selected {len(markers)} unique neutral positions -> {args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sheet", default="Sheet1")
    parser.add_argument("--header-row", default="auto", help="auto-detect coordinate headers, or a positive 1-based row number")
    parser.add_argument("--input-format", choices=("auto", "table", "positions", "bed"), default="auto")
    parser.add_argument("--chrom-column", default="1", help="1-based column number or header name")
    parser.add_argument("--pos-column", default="2", help="1-based column number or header name")
    parser.add_argument("--marker-label-column", "--label-column", dest="marker_label_column", default="TYPE",
                        help="annotation header name (case-insensitive), or 1-based column number")
    parser.add_argument("--neutral-only", action="store_true",
                        help="table already contains only neutral markers; skip annotation-column filtering")
    parser.add_argument("--neutral-marker-label", "--neutral-value", dest="neutral_marker_label", default="NEUTRAL",
                        help="keep rows matching this text in --marker-label-column; ignores case and surrounding whitespace")
    normalize(parser.parse_args())


if __name__ == "__main__":
    main()
