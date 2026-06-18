from operator import gt
import os
import gzip

import pysam


def read_fai(fai):
    fai_dict = dict()
    with open(fai, "r") as f:
        for line in f:
            if line[0] != "#":
                line = line.rstrip("\n").split("\t")
                fai_dict[line[0]] = [int(i) for i in line[1:]]
    return fai_dict


def pos_convert(p, sinf, bpl, cpl):
    return sinf + (p // bpl) * cpl + p % bpl


def get_ref(ref, fai_dict, chrom, start, end=None):
    """
    0-based coordinates, [start, end)
    ref: file handle/object of the reference
    fai_dict: dict, key: chrom, value: [chrom_len, start_in_file, base_per_line, chr_per_line]
    """

    if end is None:
        end = fai_dict[chrom][0]

    if end > fai_dict[chrom][0]:
        raise ValueError("end point exceeds the chromosome length")
    if start > fai_dict[chrom][0]:
        raise ValueError("start point exceeds the chromosome length")
    if start > end:
        raise ValueError("start point larger than end point")

    s = pos_convert(start, *fai_dict[chrom][1:])
    e = pos_convert(end, *fai_dict[chrom][1:])

    ref.seek(s)
    seq = ref.read(e-s).replace("\n", "").upper()

    return seq


def open_variant_file(vcf, mode="r", header=None):
    if header is None:
        return pysam.VariantFile(vcf, mode)
    return pysam.VariantFile(vcf, mode, header=header)


def open_text_auto(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# For now we assume the input VCF is biallelic and with sequences in the ALT field.
def infer_sv_type(ref, alts, info):
    sv_type = info.get("SVTYPE")
    if sv_type:
        return sv_type.upper()

    symbolic_types = {
        allele[1:-1].upper()
        for allele in alts
        if allele.startswith("<") and allele.endswith(">")
    }
    if len(symbolic_types) == 1:
        return next(iter(symbolic_types))

    if len(alts) == 1:
        if len(alts[0]) > len(ref):
            return "INS"
        if len(ref) > len(alts[0]):
            return "DEL"

    return None

def format_gt(sample_data):
    gt = sample_data.get("GT", (".", "."))
    alleles = [str(allele) if allele is not None else "." for allele in gt]
    separator = "|" if sample_data.phased else "/"
    formatted_gt = []
    for i, allele in enumerate(alleles):
        if i > 0:
            formatted_gt.append(separator)
        formatted_gt.append(allele)
    return formatted_gt

def parse_pysam_vcf_record(vcf_record):
    alts = list(vcf_record.alts or [])
    info = dict(vcf_record.info)
    samples = list(vcf_record.samples)
    gt = None
    if samples:
        gt = format_gt(vcf_record.samples[samples[0]])

    return {
        # "pysam_record": vcf_record,
        "chrom": vcf_record.chrom,
        "id": vcf_record.id or ".",
        "ref": vcf_record.ref,
        "alts": list(vcf_record.alts or []),
        # "filter": list(vcf_record.filter.keys()),
        # "info": info,
        "sv_type": infer_sv_type(vcf_record.ref, alts, info),
        "start": vcf_record.start,
        "end": vcf_record.start + len(vcf_record.ref),
        "alleles": vcf_record.alleles,
        # "samples": samples,
        "gt": gt,
    }


def SimpleVCFParser(vcf):
    if isinstance(vcf, (str, bytes, os.PathLike)):
        with open_variant_file(vcf, "r") as vcf_f:
            for vcf_record in vcf_f:
                yield parse_pysam_vcf_record(vcf_record)
    else:
        for vcf_record in vcf:
            yield parse_pysam_vcf_record(vcf_record)
