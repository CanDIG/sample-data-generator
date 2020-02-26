#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate random vcf files based on resampled dbSNP data

Usage:
  generate_vcf.py <number_of_files> <output_folder> <filename_start>

Options:
  -h --help          Show this screen
  <number_of_files>  Number of samples to generate
  <output_folder>    Path information where the generated files will be put
  <filename_start>   All filename will start with this phrase and continued by
                     a number.

Examples:
  generate_vcf.py 5000 tmp/variants/ patient_
"""

import os
import numpy as np
from docopt import docopt

# Path of the filtered dbSNP file - contains 1,000,000 variants
DB_FILE = r'./dbSNP_sample_1m.vcf'

SIMPLE_VCF_TEMPLATE = """##fileformat=VCFv4.0
##fileDate=20140813
##source=dbSNP
##dbSNP_BUILD_ID=141
##reference=GRCh37.p13
##phasing=partial
##INFO=<ID=RS,Number=1,Type=Integer,Description="dbSNP ID (i.e. rs number)">
##INFO=<ID=RSPOS,Number=1,Type=Integer,Description="Chr position reported in dbSNP">
##INFO=<ID=RV,Number=0,Type=Flag,Description="RS orientation is reversed">
##INFO=<ID=VP,Number=1,Type=String,Description="Variation Property.  Documentation is at ftp://ftp.ncbi.nlm.nih.gov/snp/specs/dbSNP_BitField_latest.pdf">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
##INFO=<ID=dbSNPBuildID,Number=1,Type=Integer,Description="First dbSNP Build for RS">
##INFO=<ID=SAO,Number=1,Type=Integer,Description="Variant Allele Origin: 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both">
##INFO=<ID=SSR,Number=1,Type=Integer,Description="Variant Suspect Reason Codes (may be more than one value added together) 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other">
##INFO=<ID=WGT,Number=1,Type=Integer,Description="Weight, 00 - unmapped, 1 - weight 1, 2 - weight 2, 3 - weight 3 or more">
##INFO=<ID=VC,Number=1,Type=String,Description="Variation Class">
##INFO=<ID=PM,Number=0,Type=Flag,Description="Variant is Precious(Clinical,Pubmed Cited)">
##INFO=<ID=TPA,Number=0,Type=Flag,Description="Provisional Third Party Annotation(TPA) (currently rs from PHARMGKB who will give phenotype data)">
##INFO=<ID=PMC,Number=0,Type=Flag,Description="Links exist to PubMed Central article">
##INFO=<ID=S3D,Number=0,Type=Flag,Description="Has 3D structure - SNP3D table">
##INFO=<ID=SLO,Number=0,Type=Flag,Description="Has SubmitterLinkOut - From SNP->SubSNP->Batch.link_out">
##INFO=<ID=NSF,Number=0,Type=Flag,Description="Has non-synonymous frameshift A coding region variation where one allele in the set changes all downstream amino acids. FxnClass = 44">
##INFO=<ID=NSM,Number=0,Type=Flag,Description="Has non-synonymous missense A coding region variation where one allele in the set changes protein peptide. FxnClass = 42">
##INFO=<ID=NSN,Number=0,Type=Flag,Description="Has non-synonymous nonsense A coding region variation where one allele in the set changes to STOP codon (TER). FxnClass = 41">
##INFO=<ID=REF,Number=0,Type=Flag,Description="Has reference A coding region variation where one allele in the set is identical to the reference sequence. FxnCode = 8">
##INFO=<ID=SYN,Number=0,Type=Flag,Description="Has synonymous A coding region variation where one allele in the set does not change the encoded amino acid. FxnCode = 3">
##INFO=<ID=U3,Number=0,Type=Flag,Description="In 3' UTR Location is in an untranslated region (UTR). FxnCode = 53">
##INFO=<ID=U5,Number=0,Type=Flag,Description="In 5' UTR Location is in an untranslated region (UTR). FxnCode = 55">
##INFO=<ID=ASS,Number=0,Type=Flag,Description="In acceptor splice site FxnCode = 73">
##INFO=<ID=DSS,Number=0,Type=Flag,Description="In donor splice-site FxnCode = 75">
##INFO=<ID=INT,Number=0,Type=Flag,Description="In Intron FxnCode = 6">
##INFO=<ID=R3,Number=0,Type=Flag,Description="In 3' gene region FxnCode = 13">
##INFO=<ID=R5,Number=0,Type=Flag,Description="In 5' gene region FxnCode = 15">
##INFO=<ID=OTH,Number=0,Type=Flag,Description="Has other variant with exactly the same set of mapped positions on NCBI refernce assembly.">
##INFO=<ID=CFL,Number=0,Type=Flag,Description="Has Assembly conflict. This is for weight 1 and 2 variant that maps to different chromosomes on different assemblies.">
##INFO=<ID=ASP,Number=0,Type=Flag,Description="Is Assembly specific. This is set if the variant only maps to one assembly">
##INFO=<ID=MUT,Number=0,Type=Flag,Description="Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources">
##INFO=<ID=VLD,Number=0,Type=Flag,Description="Is Validated.  This bit is set if the variant has 2+ minor allele count based on frequency or genotype data.">
##INFO=<ID=G5A,Number=0,Type=Flag,Description=">5% minor allele frequency in each and all populations">
##INFO=<ID=G5,Number=0,Type=Flag,Description=">5% minor allele frequency in 1+ populations">
##INFO=<ID=HD,Number=0,Type=Flag,Description="Marker is on high density genotyping kit (50K density or greater).  The variant may have phenotype associations present in dbGaP.">
##INFO=<ID=GNO,Number=0,Type=Flag,Description="Genotypes available. The variant has individual genotype (in SubInd table).">
##INFO=<ID=KGValidated,Number=0,Type=Flag,Description="1000 Genome validated">
##INFO=<ID=KGPhase1,Number=0,Type=Flag,Description="1000 Genome phase 1 (incl. June Interim phase 1)">
##INFO=<ID=KGPilot123,Number=0,Type=Flag,Description="1000 Genome discovery all pilots 2010(1,2,3)">
##INFO=<ID=KGPROD,Number=0,Type=Flag,Description="Has 1000 Genome submission">
##INFO=<ID=OTHERKG,Number=0,Type=Flag,Description="non-1000 Genome submission">
##INFO=<ID=PH3,Number=0,Type=Flag,Description="HAP_MAP Phase 3 genotyped: filtered, non-redundant">
##INFO=<ID=CDA,Number=0,Type=Flag,Description="Variation is interrogated in a clinical diagnostic assay">
##INFO=<ID=LSD,Number=0,Type=Flag,Description="Submitted from a locus-specific database">
##INFO=<ID=MTP,Number=0,Type=Flag,Description="Microattribution/third-party annotation(TPA:GWAS,PAGE)">
##INFO=<ID=OM,Number=0,Type=Flag,Description="Has OMIM/OMIA">
##INFO=<ID=NOC,Number=0,Type=Flag,Description="Contig allele not present in variant allele list. The reference sequence allele at the mapped position is not present in the variant allele list, adjusted for orientation.">
##INFO=<ID=WTD,Number=0,Type=Flag,Description="Is Withdrawn by submitter If one member ss is withdrawn by submitter, then this bit is set.  If all member ss' are withdrawn, then the rs is deleted to SNPHistory">
##INFO=<ID=NOV,Number=0,Type=Flag,Description="Rs cluster has non-overlapping allele sets. True when rs set has more than 2 alleles from different submissions and these sets share no alleles in common.">
##FILTER=<ID=NC,Description="Inconsistent Genotype Submission For At Least One Sample">
##INFO=<ID=CAF,Number=.,Type=String,Description="An ordered, comma delimited list of allele frequencies based on 1000Genomes, starting with the reference allele followed by alternate alleles as ordered in the ALT column. Where a 1000Genomes alternate allele is not in the dbSNPs alternate allele set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previuosly reported in VCF as the GMAF.  This is the GMAF reported on the RefSNP and EntrezSNP pages and VariationReporter">
##INFO=<ID=COMMON,Number=1,Type=Integer,Description="RS is a common SNP.  A common SNP is one that has at least one 1000Genomes population with a minor allele of frequency >= 1% and for which 2 or more founders contribute to that minor allele frequency.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
{variants}"""



def generate_length_vector(number_of_sample=3000, mu=0, sigma=0.6):
    """
    Number of somatic variants called by Strelka is on the order of magnitue:
    5k - 200k. This function samples a lognormal distribution to figure out the
    number of variants to generate for each sample.
    """
    minimum_number_of_variants = 5000
    scale = 15000

    distribution = np.random.lognormal(mu, sigma, number_of_sample)

    return (distribution * scale + minimum_number_of_variants).astype(np.uint32)


def sample_distribution(start, stop, length):
    """
    Generate unique indexes for sampling a dataset.

    np.random.choice is super slow for these many data, in contrary randint
    can generate duplicates.
    """
    indexes = np.unique(
        np.sort(
            np.random.randint(start, stop, int(length*1.1))
            )
        )[:length]

    return indexes


def sample_dbsnp(number_of_variants=1000000):
    """
    dbSNP contains almost 29 million variants. To generate variants that samples
    can share and all variants will be unique this function selects a subset of
    all dbSNP variants from which each sample will select.
    """
    # Generate variant lines to take
    num = sample_distribution(150, 28929335, number_of_variants)

    # NOTE this file is over 6 GB!
    lines = read_variant_db_file('/home/zbozoky/prj/vcf/dbSNP/common_all.vcf')

    lines = np.array(lines)[list(num)]

    with open(DB_FILE, 'w') as out:
        out.write(''.join(lines))


def read_variant_db_file(filepath):
    """
    Load dbSNP datafile
    """
    with open(filepath, 'r') as db_file:
        lines = db_file.readlines()

    return np.array(lines)


def generate_vcf(filename, variants):
    """
    It creates a single VCF file using vcf template and the variants extracted
    from dbSNP
    """
    with open(filename, 'w') as vcf_file:
        vcf_file.write(SIMPLE_VCF_TEMPLATE.format(
            variants=''.join(variants)
            ))


def generate_variant_index(vcf_filename):
    """
    Runs the following commands:

    1) SORTING
    /home/rcorbett/bin/vcftools_0.1.12a/bin/vcf-sort {VCF_filename} > {VCF_sorted}

    2) COMPRESSION
    /home/rcorbett/bin/tabix-0.2.6/bgzip -c {VCF_sorted} > {VCF.gz}

    3) INDEX FILE GENERATION
    /home/rcorbett/bin/tabix-0.2.6/tabix -f {VCF.gz}
    """
    sorter = '/home/rcorbett/bin/vcftools_0.1.12a/bin/vcf-sort'
    compressor = '/home/rcorbett/bin/tabix-0.2.6/bgzip'
    indexer = '/home/rcorbett/bin/tabix-0.2.6/tabix'


    vcf_sorted = '{0}_sorted'.format(vcf_filename)
    vcf_gzipped = '{0}.gz'.format(vcf_filename)
    vcf_index = '{0}.gz.tbi'.format(vcf_filename)


    os.system('{sorter} {VCF_filename} > {VCF_sorted}'.format(
        sorter=sorter,
        VCF_filename=vcf_filename,
        VCF_sorted=vcf_sorted,
        ))

    os.system('{compressor} -c {VCF_sorted} > {VCF_gz}'.format(
        compressor=compressor,
        VCF_sorted=vcf_sorted,
        VCF_gz=vcf_gzipped,
        ))

    os.system('{indexer} -f {VCF_gz}'.format(
        indexer=indexer,
        VCF_gz=vcf_gzipped,
        ))

    os.system('rm -rf {VCF_file} {VCF_sorted}'.format(
        VCF_file=vcf_filename,
        VCF_sorted=vcf_sorted,
        ))

    return vcf_gzipped, vcf_index



def main():
    """
    
    """
    # Deal with the arguments
    arguments = docopt(__doc__)
    number_of_files_to_generate = int(arguments['<number_of_files>'])
    output_folder = arguments['<output_folder>']
    filename_start = arguments['<filename_start>']

    # Load dbSNP information
    input_db = read_variant_db_file(DB_FILE)

    # Figure out how many variants to generate per vcf file
    vcf_lengths = generate_length_vector(number_of_files_to_generate)

    # Generate a vcf file for each sample
    for index, length in enumerate(vcf_lengths):

        # Sample the db file
        variant_indexes_to_include = sample_distribution(
            start=0,
            stop=len(input_db),
            length=length
            )
        variant_to_include = input_db[list(variant_indexes_to_include)]

        # Filename definition
        vcf_filename = '{folder}/{filename_start}{index}.vcf'.format(
            folder=output_folder,
            filename_start=filename_start,
            index=str(index+1).rjust(5, '0'),
            )

        # Save vcf file
        generate_vcf(
            filename=vcf_filename,
            variants=variant_to_include,
            )

        # Compress it and generate and index file for load
        generate_variant_index(vcf_filename)


if __name__ == '__main__':
    main()
