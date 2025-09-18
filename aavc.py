from datetime import datetime
from psycopg2.extensions import register_adapter, AsIs
from sqlalchemy import create_engine, text, pool
from typing import Any, Dict, Optional, Union, Tuple, List
import argparse
import json
import numpy
import os
import math
import pandas as pd
import re
import requests
import subprocess
import time
import traceback
import urllib3

# suppress pandas warnings about chained assignment operations and urllib3
pd.options.mode.chained_assignment = None
pd.set_option('mode.sim_interactive', True)
pd.set_option('future.no_silent_downcasting', True)
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# database connection string for PostgreSQL
with open("db_config.txt", "r") as f:
    db_url = f.read().strip()

# open a connection to the database
engine = create_engine(db_url, pool_size=200, max_overflow=0)
connection = engine.connect()

# fix issues with inserting NumPy data into PostgreSQL
register_adapter(numpy.float64, lambda x: AsIs(x))
register_adapter(numpy.int64, lambda x: AsIs(x))


class Query:
    """
    Class to process genetic variant queries, normalize input,
    and fetch annotations from VEP and gnomAD.
    """

    def __init__(self,
                 variants: List[str],
                 predictor: Any,
                 input_type: Optional[str] = None,
                 debug: bool = False,
                 cache: bool = False) -> None:
        """
        Initialize a Query object.

        Args:
            variants (List[str]): A list of variant strings (e.g., "c.123A>T").
            predictor (Any): Model or function to predict variant effects.
            input_type (Optional[str], optional): Input format type (e.g., "nt", "aa"). Defaults to None.
            debug (bool, optional): Enable debug mode. Defaults to False.
            cache (bool, optional): Enable caching of results. Defaults to False.
        """
        self.debug = debug
        self.cache = cache
        self.predictor = predictor
        self.variants = variants
        self.input_type = input_type or self.infer_input_type()

        self.normalize_input_string()
        self.eff_data: Dict[str, Any] = self.VEP()
        self.freq_data: Dict[str, Any] = self.gnomAD()

    def infer_input_type(self):

        test_input = self.variants[0]

        if re.match(r'^(chr)?([1-9]|1[0-9]|2[0-2]|(?:[xyXYmM]|MT))[-:\s][0-9]+[-:\s][-.A-Za-z]+[-:\s][-.A-Za-z]+$',
                    test_input):
            input_type = "vcf_string"
        elif "ENS" in test_input:
            input_type = "hgvs_ensembl"
        elif "NM_" in test_input:
            input_type = "hgvs_refseq"
        elif ":c" in test_input:
            input_type = "hgvs_symbol_nt"
        elif ":p" in test_input:
            input_type = "hgvs_symbol_aa"
        else:
            raise Exception(f"Invalid input: {test_input}! Example inputs: 1 54999187 A T, 1-54999187-A-T, 1:54999187-A-T, ENST00000651561:c.1A>T, NM_057176:c.1A>T, BSND:c.1A>T")

        if self.debug:
            print(input_type)

        return input_type

    def normalize_input_string(self):

        # correction for 0-based indels
        for i in range(len(self.variants)):

            if self.input_type == "vcf_string":

                # fix indels
                if "---" in self.variants[i]:
                    self.variants[i] = self.variants[i].replace('---', '-.-')
                elif "--" in self.variants[i]:
                    self.variants[i] = self.variants[i].replace('--', '-.')

                # rewrite vcf string split by dash
                tokens = re.split(r'[-:\s>\t]', self.variants[i])
                self.variants[i] = '-'.join(tokens)

            elif self.input_type in ["hgvs_ensembl", "hgvs_refseq"]:

                # omit version number, vep cannot read them
                self.variants[i] = re.sub(r'\.\d+:', ':', self.variants[i])

    @staticmethod
    def convert_nt_notation(notation):
        if "+" in notation or "-" in notation or "_" in notation:
            return notation

        else:
            pattern = r'(\d+)([ACGTacgt])>([ACGTacgt])'
            converted = re.sub(pattern, r'\2\1\3', notation)
            return converted

    @staticmethod
    def convert_aa_notation(notation):

        amino_acid_mapping = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
            'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G',
            'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
            'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
            'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
            'Ter': 'X'
        }
        # iterate through the input string and convert three-letter codes to single-letter codes
        result = ''
        i = 0
        while i < len(notation):
            found = False
            for code, letter in amino_acid_mapping.items():
                if notation[i:i + 3] == code:
                    result += letter
                    i += 3
                    found = True
                    break
            if not found:
                result += notation[i]
                i += 1

        return result

    def VEP(self):

        query_list = []
        var_responses = {}

        for variant in self.variants:
            if self.input_type == "vcf_string":

                # if indel, recode vcf-string
                if "." in variant:

                    split = variant.split("-")
                    chr = split[0]
                    pos = int(split[1])
                    ref = split[2]
                    alt = split[3]

                    # convert vcf-string to vep input
                    if ref == ".":
                        ref = "-"
                    elif alt == ".":
                        alt = "-"

                    # calculate start and end positions
                    if ref == "-":
                        start_pos = pos
                        end_pos = pos - 1
                    elif len(ref) == 1:
                        start_pos = pos
                        end_pos = pos
                    else:
                        start_pos = pos
                        end_pos = pos + len(ref) - 1

                    # add the decoded indel into the list
                    query_list.append(f"{chr} {start_pos} {end_pos} {ref}/{alt} 1")
                else:
                    split = variant.split("-")
                    chr = split[0]
                    pos = int(split[1])
                    ref = split[2]
                    alt = split[3]

                    # add the decoded snp into the list
                    query_list.append(f"{chr} {pos} . {ref} {alt} . . .")
            else:
                query_list.append(variant)

        vep_query = '{{ "variants" : {0} }}'.format(json.dumps(query_list))

        print(query_list)

        server = "https://rest.ensembl.org"
        ext = "/vep/homo_sapiens/region"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}

        attempts = 0
        # retry if connection lost or no response received
        while attempts < 16:

            print("retrieving via REST API")
            # select prioritization strategy


            if self.predictor == "revel":
                pred_field = "REVEL_score"
            else:
                pred_field = "BayesDel_noAF_score"

            response = requests.post(server + ext, headers=headers, data=vep_query,
                                     params={"pick": 1, "hgvs": 1,
                                             "pick_order": "rank,mane_select,canonical,tsl,biotype,ccds,length",
                                             #"pick_order": "mane_select,canonical,tsl,biotype,ccds,rank,length",
                                             "dbNSFP": f"phyloP100way_vertebrate,{pred_field}",
                                             "SpliceAI": 1, "vcf_string": 1, "canonical": 1, "numbers": 1})

            if response.ok:
                decoded = response.json()

                break

            elif response.status_code == 400:
                raise Exception("Invalid input!")
            else:
                attempts += 1
                print(response)
                print(f"{response} Could not establish connection with Ensembl. Retrying...")
                time.sleep(attempts * 1.5)

        else:
            raise Exception("Could not access to Ensembl. Aborted.")

        i = 0
        for var in decoded:

            if self.debug:
                print(var)

            try:

                # select the prioritized variant effect
                j = 0


                strand = var["transcript_consequences"][j]["strand"]
                gene_id = var["transcript_consequences"][j]["gene_id"]
                if "gene_symbol" in var["transcript_consequences"][j]:
                    gene = var["transcript_consequences"][j]["gene_symbol"]
                else:
                    gene = gene_id
                tcpt_id = var["transcript_consequences"][j]["transcript_id"]
                major_conseq = var["most_severe_consequence"]
                all_conseqs = var["transcript_consequences"][j]["consequence_terms"]

                if self.predictor == "revel":

                    if "revel_score" in var["transcript_consequences"][j]:
                        revel_raw = var["transcript_consequences"][j]["revel_score"]

                        if isinstance(revel_raw, str):
                            all_scores = revel_raw.split(',')
                            revel = None

                            for score in all_scores:
                                if score != ".":
                                    revel = float(score)
                                    break
                        else:
                            revel = float(revel_raw)
                    else:
                        revel = None

                    bayesdel = None

                elif self.predictor == "bayesdel":

                    if "bayesdel_noaf_score" in var["transcript_consequences"][j]:
                        bayesdel_raw = var["transcript_consequences"][j]["bayesdel_noaf_score"]

                        if isinstance(bayesdel_raw, str):
                            all_scores = bayesdel_raw.split(',')
                            bayesdel = None

                            for score in all_scores:
                                if score != ".":
                                    bayesdel = float(score)
                                    break
                        else:
                            bayesdel = float(bayesdel_raw)
                    else:
                        bayesdel = None

                    revel = None
                if "phylop100way_vertebrate" in var["transcript_consequences"][j]:
                    phyloP100way = var["transcript_consequences"][j]["phylop100way_vertebrate"]
                else:
                    phyloP100way = None

                if "hgvsc" in var["transcript_consequences"][j]:
                    nt_change = self.convert_nt_notation(
                        var["transcript_consequences"][j]["hgvsc"].split(":")[1])
                else:
                    nt_change = None

                if "hgvsp" in var["transcript_consequences"][j]:
                    aa_change = self.convert_aa_notation(
                        var["transcript_consequences"][j]["hgvsp"].split(":")[1])
                    aa_pos = int(re.match(r"p\.[A-Za-z]([0-9]+).*", aa_change).group(1))
                else:
                    aa_change = None
                    aa_pos = None

                if "spliceai" in var["transcript_consequences"][j]:

                    DS_AG = var["transcript_consequences"][j]["spliceai"]["DS_AG"]
                    DS_AL = var["transcript_consequences"][j]["spliceai"]["DS_AL"]
                    DS_DG = var["transcript_consequences"][j]["spliceai"]["DS_DG"]
                    DS_DL = var["transcript_consequences"][j]["spliceai"]["DS_DL"]
                    '''
                    DP_AG = var["transcript_consequences"][j]["spliceai"]["DP_AG"]
                    DP_AL = var["transcript_consequences"][j]["spliceai"]["DP_AL"]
                    DP_DG = var["transcript_consequences"][j]["spliceai"]["DP_DG"]
                    DP_DL = var["transcript_consequences"][j]["spliceai"]["DP_DL"]
                    '''
                    max_spl = max(DS_AG, DS_AL, DS_DG, DS_DL)
                else:
                    max_spl = None

                is_canonical = False
                if "canonical" in var["transcript_consequences"][j]:
                    if var["transcript_consequences"][j]["canonical"] == 1:
                        is_canonical = True

                # collect decoded info
                var_data = {"gene": gene, "gene_id": gene_id, "tcpt_id": tcpt_id, "strand": strand,
                            "nt_change": nt_change, "aa_change": aa_change, "aa_pos": aa_pos,
                            "major_conseq": major_conseq, "all_conseqs": all_conseqs, "revel": revel,
                            "bayesdel": bayesdel, "spliceAI": max_spl, "phyloP100way": phyloP100way,
                            "is_canonical": is_canonical}

                print(var_data)
                # if indel or non-vcf-string identifier used, adopt 1-based positions
                if "." in self.variants[i] or self.input_type != "vcf_string":
                    print(var["vcf_string"])
                    self.variants[i] = var["vcf_string"]

            except Exception as e:
                print(f"Could not decipher: {self.variants[i]}: {e}")
                var_data = {"gene": None, "gene_id": None, "tcpt_id": None, "strand": None,
                            "nt_change": None, "aa_change": None, "aa_pos": None, "major_conseq": None,
                            "all_conseqs": None, "revel": None, "bayesdel": None, "spliceAI": None,
                            "phyloP100way": None, "is_canonical": None}

            var_responses[self.variants[i]] = var_data
            i += 1

        return var_responses

    def gnomAD(self):

        gnomAD_vars = {}
        query_list = []

        i = 1
        for var_id in self.variants:
            print(f"mapping: {var_id}")
            varqry = f'''variant{i}: variant(variantId: "{var_id}", dataset: gnomad_r4) {{exome {{homozygote_count hemizygote_count filters faf95 {{popmax}}}}}}'''
            query_list.append(varqry)
            i += 1

        fmt = "{" + " ".join(query_list) + "}"

        retry_count = 0
        while retry_count < 32:

            response = requests.post("https://gnomad.broadinstitute.org/api",
                                     json={"query": fmt, "variables": {}},
                                     headers={"Content-Type": "application/json"})

            if response.status_code == 200:
                decoded = response.json()["data"]

                j = 1
                for var_id in self.variants:

                    var_data = decoded[f"variant{j}"]

                    if var_data and var_data["exome"]:

                        if var_data["exome"]["filters"]:

                            variant = {"popmax": 0, "total_hom_ac": 0}

                        else:

                            total_hom_ac = var_data["exome"]["homozygote_count"] + \
                                           var_data["exome"]["hemizygote_count"]

                            popmax = var_data["exome"]["faf95"]["popmax"]

                            if not popmax:
                                popmax = 0

                            variant = {"popmax": popmax, "total_hom_ac": total_hom_ac}
                    else:

                        variant = {"popmax": 0, "total_hom_ac": 0}

                    gnomAD_vars[var_id] = variant
                    j += 1

                print("\n")
                return gnomAD_vars
            else:
                retry_count += 1
                time.sleep(retry_count * 1.5)
                print(f"Retry #{retry_count} - Status Code: {response.status_code}")


class Variant:
    """
    Class representing a single genetic variant with effect, frequency,
    transcript, and disease-specific annotations.
    """

    def __init__(self,
                 var_id: str,
                 eff_data: Dict[str, Any],
                 freq_data: Dict[str, Any],
                 predictor: Any) -> None:
        """
        Initialize a Variant object.

        Args:
            var_id (str): Variant identifier in the format "chr-pos-ref-alt".
            eff_data (Dict[str, Any]): Variant effect annotation (e.g., VEP output).
            freq_data (Dict[str, Any]): Variant frequency information (e.g., gnomAD).
            predictor (Any): Predictor object or model for functional effect.
        """
        # query variables
        self.var_id: str = var_id
        self.eff_data: Dict[str, Any] = eff_data
        self.freq_data: Dict[str, Any] = freq_data
        self.predictor: Any = predictor

        # basic variant info
        split = var_id.split("-")
        self.chr: str = split[0]
        self.pos: int = int(split[1])
        self.ref: str = split[2]
        self.alt: str = split[3]
        self.end_pos: int = self.pos + len(self.ref) - 1 if len(self.ref) > 1 else self.pos

        # genomic context
        self.gene: str = eff_data["gene"]
        self.gene_id: str = eff_data["gene_id"]
        self.tcpt_id: str = eff_data["tcpt_id"]
        self.strand: str = eff_data["strand"]
        self.nt_change: Optional[str] = eff_data.get("nt_change")

        # nucleotide position in transcript
        self.nt_pos: Optional[int] = None
        if self.nt_change and not ("+" in self.nt_change or "-" in self.nt_change):
            self.nt_pos = int(re.search(r'\d+', self.nt_change).group())

        # effect info
        self.aa_change: Optional[str] = eff_data.get("aa_change")
        self.aa_pos: Optional[int] = eff_data.get("aa_pos")
        self.REVEL: Optional[float] = eff_data.get("revel")
        self.BayesDel: Optional[float] = eff_data.get("bayesdel")
        self.spliceAI: Optional[float] = eff_data.get("spliceAI")
        self.phyloP100way: Optional[float] = eff_data.get("phyloP100way")
        self.major_conseq: str = eff_data.get("major_conseq", "")
        self.all_conseqs: List[str] = eff_data.get("all_conseqs", [])
        self.domain: Optional[str] = self.check_domains()

        # frequency info
        self.filt_af: float = freq_data.get("popmax", 0)
        self.hom_ac: int = freq_data.get("total_hom_ac", 0)

        # gene-specific metrics
        (self.pLI, self.LOEUF, self.mis_z, self.HI_score, self.mis_perc,
         self.lof_perc, self.BS1_cutoff, self.BS2_cutoff) = self.get_gene_metrics()

        # disease associations
        self.LOF_assoc, self.MOI = self.get_disease_assoc()
        self.is_in_LOF_gene: bool = self.is_in_LOF_gene()

        # transcript info
        tcpt = self.fetch_tcpt()
        self.prot_length: int = tcpt["prot_length"]
        self.is_canonical: bool = tcpt["is_canonical"]
        self.alt_start: int = tcpt["alt_start"]
        self.alt_stop: int = tcpt["alt_stop"]
        self.tln_start: int = tcpt["tln_start"]

        # splicing
        self.is_splicing: bool = "splic" in self.major_conseq or any("splic" in c for c in self.all_conseqs)
        self.splice_type: Optional[str] = None

        # exon positions
        self.exon: Union[pd.DataFrame, None] = self.find_exon()
        if isinstance(self.exon, pd.DataFrame) and not self.exon.empty:
            self.tcpt_start = self.exon["tcpt_start"].squeeze()
            self.tcpt_end = self.exon["tcpt_end"].squeeze()
            self.exon_start = self.exon["exon_start"].squeeze()
            self.exon_end = self.exon["exon_end"].squeeze()
            self.coding_end = self.exon["coding_end"].squeeze()
            self.highest_freq = float(self.exon["highest_freq"].squeeze() or 0)
        else:
            self.tcpt_start = self.tcpt_end = self.exon_start = self.exon_end = self.coding_end = None
            self.highest_freq = 0

        # NMD prediction
        self.undergoes_NMD: bool = False if self.is_at_canonical_splice_site() else self.check_NMD()

        # ClinVar info
        self.is_denovo: int = 0
        self.is_intrans: int = 0
        self.clinvar: Any = self.check_clinvar_status()

        # ACMG classification placeholders
        self.ACMG_criteria: Optional[Dict[str, Any]] = None
        self.ACMG_classification: Optional[str] = None
        self.ACMG_score: Optional[float] = None
        self.ACMG_flags: Optional[List[str]] = None
        self.post_prob: float = 0.0

    def fetch_tcpt(self):

        query_tcpt = text(''' SELECT * FROM exons WHERE enst_id = :tcpt_id ''')
        tcpt_row = pd.read_sql_query(query_tcpt, connection, params={"tcpt_id": self.tcpt_id})

        if len(tcpt_row) > 0:

            strand = tcpt_row["strand"][0].squeeze()
            is_canonical = tcpt_row["is_canonical"][0].squeeze()
            alt_start = tcpt_row["alt_start"][0].squeeze()
            alt_stop = tcpt_row["alt_stop"][0].squeeze()
            prot_length = tcpt_row["prot_length"][0]
            tln_start = tcpt_row["tss"][0]

            if prot_length:
                prot_length = prot_length.squeeze()
            else:
                prot_length = None

        else:
            # it could continue with canonical transcript belonging to the same gene
            # on a second thought, I concluded that it would to lead to errors
            # considering all tcpt matching problems we see are due to version differences
            # the program will raise an exception for the variant and skip it

            if len(tcpt_row) > 0:
                raise Exception("Transcript/gene could not be deciphered")
            else:
                raise Exception("Transcript/gene could not be matched")

        return {"strand": strand, "prot_length": prot_length, "tln_start": tln_start,
                "is_canonical": is_canonical, "alt_start": alt_start, "alt_stop": alt_stop}

    def check_clinvar_status(self):

        q = text(''' SELECT * FROM clinvar WHERE chr = :chr AND pos = :pos AND ref = :ref AND alt = :alt ''')
        r = pd.read_sql_query(q, connection,
                              params={"chr": self.chr, "pos": self.pos, "ref": self.ref, "alt": self.alt})

        if len(r) == 1:

            self.is_intrans = r["is_intrans"].squeeze()
            self.is_denovo = r["is_denovo"].squeeze()

        elif len(r) > 1:

            if len(r.loc[r["is_intrans"] == 1, 'is_intrans'].squeeze()):
                self.is_intrans = 1
            else:
                self.is_intrans = 0

            if len(r.loc[r["is_denovo"] == 1, 'is_denovo'].squeeze()):
                self.is_denovo = 1
            else:
                self.is_denovo = 0

        return r

    def clinvar_within(self, method="exon", query_start=None, query_end=None):

        # if method is gene or splicing, function returns its own result
        # otherwise, based on method, query_start and end positions are set

        if method == "gene":

            q_clv = text(''' SELECT
                                    COUNT(CASE WHEN cls = 'P' OR "cls" = 'LP'
                                        THEN 1 END) as pathogenic,
                                    COUNT(CASE WHEN cls = 'B' OR "cls" = 'LB'
                                        THEN 1 END) as benign
                                    FROM clinvar WHERE gene LIKE :gene AND is_missense = 1 ''')

            r_clv = connection.execute(q_clv, {"gene": f"%{self.gene}%"}).fetchone()

            return r_clv

        elif method == "splice_site":

            # check if variant is located canonical splice sites
            # whether it is at donor or acceptor, empirically inferred
            if self.is_at_canonical_splice_site:
                if self.strand == 1:

                    q1 = text('''SELECT * FROM exons WHERE enst_id = :tcpt_id AND :pos - exon_end IN (1, 2)''')
                    q2 = text('''SELECT * FROM exons WHERE enst_id = :tcpt_id AND exon_start - :pos IN (1, 2)''')

                else:

                    q1 = text('''SELECT * FROM exons WHERE enst_id = :tcpt_id AND exon_end - :pos IN (1, 2)''')
                    q2 = text('''SELECT * FROM exons WHERE enst_id = :tcpt_id AND :pos - exon_start IN (1, 2)''')

            else:
                if self.strand == 1:

                    q1 = text('''SELECT * FROM exons WHERE enst_id = :tcpt_id AND 
                                ((exon_end - :pos BETWEEN 0 AND 2) OR (:pos - exon_end BETWEEN 3 AND 6))''')

                    q2 = text('''SELECT * FROM exons WHERE enst_id = :tcpt_id AND 
                                (exon_start = :pos OR exon_start - :pos BETWEEN 3 AND 20)''')
                else:

                    q1 = text('''SELECT * FROM exons WHERE enst_id = :tcpt_id AND 
                                    ((:pos - exon_end BETWEEN 0 AND 2) OR (exon_end - :pos BETWEEN 3 AND 6))''')

                    q2 = text('''SELECT * FROM exons WHERE enst_id = :tcpt_id AND 
                                    (exon_start = :pos OR :pos - exon_start BETWEEN 3 AND 20)''')

            donor = pd.read_sql_query(q1, connection, params={"pos": self.pos, "tcpt_id": self.tcpt_id})
            acceptor = pd.read_sql_query(q2, connection, params={"pos": self.pos, "tcpt_id": self.tcpt_id})

            if len(donor) > 0 and len(acceptor) < 1:
                exon = donor
                type = "donor"
            elif len(acceptor) > 0 and len(donor) < 1:
                exon = acceptor
                type = "acceptor"
            else:
                type = None
                exon = None

            if type == "donor":

                if self.strand == 1:

                    query_start = exon["exon_end"].squeeze() - 3
                    query_end = exon["exon_end"].squeeze() + 6

                    dint_start = exon["exon_end"].squeeze() + 1
                    dint_end = exon["exon_end"].squeeze() + 2

                else:

                    query_start = exon["exon_end"].squeeze() + 3
                    query_end = exon["exon_end"].squeeze() - 6

                    dint_start = exon["exon_end"].squeeze() - 1
                    dint_end = exon["exon_end"].squeeze() - 2

            elif type == "acceptor":

                if self.strand == 1:

                    query_start = exon["exon_start"].squeeze() - 20
                    query_end = exon["exon_start"].squeeze()

                    dint_start = exon["exon_start"].squeeze() - 2
                    dint_end = exon["exon_start"].squeeze() - 1

                else:

                    query_start = exon["exon_start"].squeeze() + 20
                    query_end = exon["exon_start"].squeeze()

                    dint_start = exon["exon_start"].squeeze() + 2
                    dint_end = exon["exon_start"].squeeze() + 1

            else:
                return None

            ran_vals = [query_start, query_end]
            query_start = min(ran_vals)
            query_end = max(ran_vals)

            q_tcpt = text(''' SELECT * FROM clinvar WHERE enst_id LIKE :tcpt_id ''')
            same_tcpt = pd.read_sql_query(q_tcpt, connection, params={"tcpt_id": f"%{self.tcpt_id}%"})

            same_nt = same_tcpt[(same_tcpt["pos"] == self.pos)]

            p_same_nt = len(same_nt[same_nt["cls"] == "P"])
            lp_same_nt = len(same_nt[same_nt["cls"] == "LP"])

            p_same_region = len(same_tcpt[(same_tcpt["pos"] <= query_start) & (same_tcpt["pos"] >= query_end)] &
                                same_tcpt[same_tcpt["cls"] == "P"])
            lp_same_region = len(same_tcpt[(same_tcpt["pos"] <= query_start) & (same_tcpt["pos"] >= query_end)] &
                                 same_tcpt[same_tcpt["cls"] == "LP"])

            p_same_dint = len(same_tcpt[(same_tcpt["pos"] <= dint_start) &
                                        (same_tcpt["pos"] >= dint_end) &
                                        (same_tcpt["cls"] == "P")])

            return [p_same_nt, lp_same_nt, p_same_region, lp_same_region, p_same_dint]

        elif method == "alternative_start":

            query_start = self.tln_start  #self.tcpt_start

            if self.strand == 1:
                query_end = self.tln_start + self.alt_start * 3 - 3
            else:
                query_end = self.tln_start - self.alt_start * 3 + 3

        elif method == "proximity":

            if self.strand == 1:
                query_start = self.pos - 30
                query_end = self.pos + 30
            else:
                query_start = self.pos + 30
                query_end = self.pos - 30

        elif method == "domain":

            query_start = self.domain["domain_start"]
            query_end = self.domain["domain_end"]

        elif method == "range":

            query_start = query_start
            query_end = query_start

        else:  # exon
            query_start = self.exon_start
            query_end = self.exon_end

        q_clv = text(''' SELECT
                              COUNT(CASE WHEN cls = 'P' OR cls = 'LP'
                                  THEN 1 END) as pathogenic,
                              COUNT(CASE WHEN cls = 'B' OR cls = 'LB'
                                  THEN 1 END) as benign
                              FROM clinvar WHERE chr = :chr AND ((pos BETWEEN :query_start AND :query_end) OR  (pos BETWEEN :query_end AND :query_start)) AND is_missense = 1''')

        r_clv = connection.execute(q_clv,
                                   {"chr": self.chr, "query_start": query_start, "query_end": query_end}).fetchone()

        return r_clv

    def check_repeat(self):

        q_rp = text(
            ''' SELECT * FROM repeats WHERE chr = :chr AND strand = :strand AND ((start_pos <= :var_start AND end_pos >= :var_start) OR (start_pos <= :var_end AND end_pos >= :var_end)) ''')

        repeats = pd.read_sql_query(q_rp, connection,
                                    params={"var_start": self.pos, "var_end": self.end_pos, "chr": self.chr,
                                            "strand": self.strand})

        within_repeat = False
        contains_pathogenic = False
        repeat_list = []

        if len(repeats) == 1:
            within_repeat = True
            if repeats["pat_vars"].squeeze() > 0:
                contains_pathogenic = True
            repeat_list.append(
                f"{repeats['rep_name'].squeeze()} ({repeats['rep_class'].squeeze()}) ({repeats['pat_vars'].squeeze()} P/LP)")
        elif len(repeats) > 1:
            within_repeat = True
            if repeats["pat_vars"].sum() > 0:
                contains_pathogenic = True
            for index, row in repeats.iterrows():
                repeat_list.append(
                    f"{row['rep_name']} ({row['rep_class']}) ({row['pat_vars']} P/LP)")

        return {"within_repeat": within_repeat, "contains_pathogenic": contains_pathogenic, "repeat_desc": repeat_list}

    def find_exon(self):

        if self.is_splicing is True:
            q_tcpt = text(''' SELECT * FROM exons WHERE enst_id = :tcpt_id ''')

            tcpt = pd.read_sql_query(q_tcpt, connection, params={"tcpt_id": self.tcpt_id})

            if self.strand == 1:

                # -3 to +6
                donor = tcpt.loc[(tcpt["exon_end"] - 2 <= self.pos) & (self.pos <= tcpt["exon_end"] + 8)]
                # +1 to -20
                acceptor = tcpt.loc[(tcpt["exon_start"] - 20 <= self.pos) & (self.pos <= tcpt["exon_start"])]
            else:
                # -3 to +6
                donor = tcpt.loc[(tcpt["exon_end"] + 2 >= self.pos) & (self.pos >= tcpt["exon_end"] - 8)]
                # +1 to -20
                acceptor = tcpt.loc[(tcpt["exon_start"] + 20 >= self.pos) & (self.pos >= tcpt["exon_start"])]
            exon = None

            if len(donor) > 0 and len(acceptor) < 1:
                exon = donor
                self.splice_type = "donor"

            elif len(acceptor) > 0 and len(donor) < 1:
                exon = acceptor
                self.splice_type = "acceptor"

            if isinstance(exon, pd.DataFrame):
                return exon

        if self.strand == 1:
            q_exon = text(" SELECT * FROM exons WHERE enst_id = :tcpt_id AND :pos BETWEEN exon_start AND exon_end ")
        else:
            q_exon = text(" SELECT * FROM exons WHERE enst_id = :tcpt_id AND :pos BETWEEN exon_end AND exon_start ")

        exon = pd.read_sql_query(q_exon, connection, params={"tcpt_id": self.tcpt_id, "pos": self.pos})

        return exon

    def is_at_canonical_splice_site(self):

        if not self.exon_start:
            return False

        else:

            rel_to_start = (self.pos - self.exon_start) * self.strand
            rel_to_end = (self.pos - self.exon_end) * self.strand

            if rel_to_start == -1 or rel_to_start == -2:

                if self.strand == -1:
                    self.splice_type = "donor"
                else:
                    self.splice_type = "acceptor"

                return True

            elif rel_to_end == 1 or rel_to_end == 2:
                if self.strand == -1:
                    self.splice_type = "acceptor"
                else:
                    self.splice_type = "donor"

                return True

            else:
                return False

    def is_blacklisted(self):

        if self.strand == 1:
            q_bl = text(''' SELECT * FROM blacklist WHERE chr = :chr AND :pos BETWEEN "start_pos" AND "end_pos"  ''')
        else:
            q_bl = text(''' SELECT * FROM blacklist WHERE chr = :chr AND :pos BETWEEN "end_pos" AND "start_pos"  ''')

        blackrow = pd.read_sql_query(q_bl, connection, params={"pos": self.pos, "chr": self.chr})

        if len(blackrow) > 0:
            return True
        else:
            return False

    def check_NMD(self, splice_trunc_start=None):

        can_undergo_NMD = False

        for conseq in self.all_conseqs:
            if conseq == "stop_gained":
                can_undergo_NMD = True
                break
            elif splice_trunc_start is not None:
                can_undergo_NMD = True
            elif conseq == "frameshift_variant":
                can_undergo_NMD = True
                break
            else:
                return False

        if can_undergo_NMD:

            exon_count = self.exon["exon_count"].squeeze()
            exon_rank = self.exon["exon_rank"].squeeze()

            if splice_trunc_start:
                rel_trunc_start = abs(splice_trunc_start - self.tln_start)
            else:
                if "_" in self.nt_change and "*" not in self.nt_change:  #not trunc_start and
                    rel_trunc_start = int(re.search(r'\d+', self.nt_change).group())
                else:
                    rel_trunc_start = abs(self.pos - self.tln_start)

            rel_pos = abs(self.coding_end - rel_trunc_start)

            if rel_trunc_start <= 100:
                return False
            elif exon_count == 1:
                return False
            elif exon_count == exon_rank:
                return False
            elif exon_count - exon_rank == 1 and 0 <= rel_pos <= 55:
                return False
            else:
                return True
        else:
            return False

    def classify(self, PS4, deactivate_PM2, deactivate_PP5_BP6):

        new_assessment = ACMG(self, predictor=self.predictor)
        new_assessment.classify(PS4, deactivate_PM2, deactivate_PP5_BP6)

    def check_domains(self, extent="variant", query_start=None, query_end=None):

        if extent == "variant":

            if not isinstance(self.aa_pos, int):
                return {"domain": None, "domain_start": None, "domain_end": None}

            q = text("SELECT * FROM domains WHERE chr = :chr AND strand = :strand AND "
                     "start_pos <= :pos AND end_pos >= :pos ORDER BY ABS(start_pos - end_pos) LIMIT 1 ")

            r = pd.read_sql_query(q, connection, params={"chr": self.chr, "strand": self.strand, "pos": self.pos})

            if len(r) > 0:
                return {"domain": r["domain_id"].squeeze(), "domain_start": r["start_pos"].squeeze(),
                        "domain_end": r["end_pos"].squeeze()}
            else:
                return None

        elif extent == "range":

            if self.strand == 1:
                q = text(" SELECT * FROM domains WHERE chr = :chr AND strand = :strand AND "
                         "(start_pos >= :query_start AND end_pos <= :query_start) OR "
                         "(start_pos >= :query_end AND end_pos <= :query_end)")
            else:
                q = text(" SELECT * FROM domains WHERE chr = :chr AND strand = :strand AND "
                         "(start_pos <= :query_start AND end_pos >= :query_start) OR "
                         "(start_pos <= :query_end AND end_pos >= :query_end)")

            r = pd.read_sql_query(q, connection, params={"chr": self.chr, "strand": self.strand,
                                                         "query_start": query_start, "query_end": query_end})

            domains = {}

            if len(r) > 0:
                for index, domain in r.iterrows():
                    domains[domain["domain_id"]] = [domain["start_pos"], domain["end_pos"]]
                return domains
            else:
                return None

        elif extent == "tcpt":

            if self.strand == 1:
                q = text(" SELECT * FROM domains WHERE chr = :chr AND strand = :strand AND "
                         "start_pos >= :tcpt_start AND end_pos >= :tcpt_start AND "
                         "start_pos <= :tcpt_end AND end_pos <= :tcpt_end")
            else:
                q = text(" SELECT * FROM domains WHERE chr = :chr AND strand = :strand AND "
                         "start_pos <= :tcpt_start AND end_pos <= :tcpt_start AND "
                         "start_pos >= :tcpt_end AND end_pos >= :tcpt_end")

            r = pd.read_sql_query(q, connection, params={"chr": self.chr, "strand": self.strand,
                                                         "tcpt_start": self.tcpt_start, "tcpt_end": self.tcpt_end})

            domains = {}

            if len(r) > 0:
                for index, domain in r.iterrows():
                    domains[domain["domain_id"]] = [domain["start_pos"], domain["end_pos"]]
                return domains
            else:
                return None


        if self.strand == 1:
            q = text(
                '''SELECT * FROM pext WHERE ensg_id = :gene_id AND chr = :chromosome AND :nt_pos BETWEEN start_pos AND end_pos''')
        else:
            q = text(
                '''SELECT * FROM pext WHERE ensg_id = :gene_id AND chr = :chromosome AND :nt_pos BETWEEN end_pos AND start_pos''')

        r = pd.read_sql_query(q, connection,
                              params={"gene_id": self.gene_id, "chromosome": self.chr, "nt_pos": self.pos})
        if len(r) > 0:
            return r["mean_pext"].squeeze()

        else:
            return None

    def check_homopolymer(self):

        if self.strand == 1:
            q = text(
                ''' SELECT * FROM homopolymers WHERE chr = :chr AND ((start_pos <= :var_start - 1 AND end_pos >= :var_start - 1) OR (start_pos <= :var_end + 1 AND end_pos >= :var_end + 1))''')
        else:
            q = text(
                ''' SELECT * FROM homopolymers WHERE chr = :chr AND ((start_pos >= :var_start + 1 AND end_pos <= :var_start + 1) OR (start_pos >= :var_end - 1 AND end_pos <= :var_end - 1))''')

        r = pd.read_sql_query(q, connection, params={"chr": self.chr, "var_start": self.pos, "var_end": self.end_pos})

        if len(r) > 0:
            return True
        else:
            return False

    def check_regional_constraint(self, query_start=None, query_end=None):

        if query_start:

            if self.strand == 1:
                q = text(" SELECT * FROM regional_constraints WHERE chr = :chr AND "
                         "start_pos >= :query_start AND end_pos >= :query_start AND "
                         "start_pos <= :query_end AND end_pos <= :query_end")
            else:
                q = text(" SELECT * FROM regional_constraints WHERE chr = :chr AND "
                         "start_pos <= :query_start AND end_pos <= :query_start AND "
                         "start_pos >= :query_end AND end_pos >= :query_end")

            r = pd.read_sql_query(q, connection,
                                  params={"chr": self.chr, "strand": self.strand, "query_start": query_start,
                                          "query_end": query_end})

        else:

            if self.strand == 1:
                q = text(
                    ''' SELECT * FROM regional_constraints WHERE chr = :chromosome AND :nt_pos BETWEEN start_pos AND end_pos''')
            else:
                q = text(
                    ''' SELECT * FROM regional_constraints WHERE chr = :chromosome AND :nt_pos BETWEEN end_pos AND start_pos''')

            r = pd.read_sql_query(q, connection, params={"chromosome": self.chr, "nt_pos": self.pos})

        if len(r) > 0:
            return True
        else:
            return False

    def get_gene_metrics(self):

        q_cons = text(''' SELECT * FROM gene_metrics WHERE enst_id = :tcpt_id ''')

        r_cons = pd.read_sql_query(q_cons, connection, params={"tcpt_id": self.tcpt_id})

        r_cons = r_cons.fillna(-1)

        if len(r_cons) > 0:
            pLI = r_cons["pli"].squeeze()
            LOEUF = r_cons["loeuf"].squeeze()
            mis_z = r_cons["mis_z"].squeeze()
            HI_score = r_cons["hi_score"].squeeze()
            mis_perc = r_cons["mis_perc"].squeeze()
            lof_perc = r_cons["lof_perc"].squeeze()
            freq = r_cons["freq_cutoff"].squeeze()
            nhom = r_cons["nhomalt"].squeeze()

            if freq > 0:

                _freq = freq

                decimals = 0
                while _freq < 1:
                    _freq *= 10
                    decimals += 1

                BS1_cutoff = math.ceil(freq * (10 ** decimals)) / (10 ** decimals)
                BS2_cutoff = nhom + 2
            else:

                BS1_cutoff = -1
                BS2_cutoff = -1

            return [pLI, LOEUF, mis_z, HI_score, mis_perc, lof_perc, BS1_cutoff, BS2_cutoff]


        else:

            return [-1, -1, -1, -1, -1, -1, -1, -1]

    def get_disease_assoc(self):

        q_assoc = text(''' SELECT * FROM disease_assoc WHERE ensg_id = :ensg_id ''')

        r_assoc = pd.read_sql_query(q_assoc, connection, params={"ensg_id": self.gene_id})

        mop_list = []
        moi_list = []

        for index, assoc in r_assoc.iterrows():

            mop = assoc["mop"]
            moi = assoc["moi"]

            if mop:
                mop_list.append(mop)

            if moi:
                moi_list.append(moi)

        mop_list = list(set(mop_list))
        moi_list = list(set(moi_list))

        if "LOF" in mop_list:
            LOF = True
        else:
            LOF = False

        if "BOTH" in moi_list or ("BIALLELIC" in moi_list and "MONOALLELIC" in moi_list):
            MOI = "BOTH"
        elif "X-LINKED" in moi_list:
            MOI = "X-LINKED"
        elif len(moi_list) == 1:
            MOI = moi_list[0]
        else:
            MOI = "UNKNOWN"

        return [LOF, MOI]

    def is_in_LOF_gene(self):

        if (((3 <= self.HI_score <= 4) or
             (self.HI_score != 40 and -1 < self.LOEUF < 0.6) or
             (self.HI_score != 40 and -1 < self.pLI >= 0.9)) or
                self.LOF_assoc is True):
            return True
        else:
            return False

    def is_domain_lost(self, query_start=None, query_end=None):

        if query_start:
            domains = self.check_domains(extent="range", query_start=query_start, query_end=query_end)
        else:
            domains = self.check_domains(extent="tcpt")
            query_start = self.pos
            query_end = self.tcpt_end

        if domains:

            for domain, domain_range in domains.items():

                domain_start = domain_range[0]
                domain_end = domain_range[1]

                if self.strand == 1:
                    if (query_start <= domain_start <= query_end) or (query_start <= domain_end <= query_end):
                        return True
                    else:
                        return False
                else:
                    if (query_start >= domain_start >= query_end) or (query_start >= domain_end >= query_end):
                        return True
                    else:
                        return False

        else:
            return False

    def get_func_class(self):

        q = text(''' SELECT * FROM functional WHERE chr = :chr AND pos = :pos AND ref = :ref AND alt = :alt''')

        r = pd.read_sql_query(q, connection,
                              params={"chr": self.chr, "pos": self.pos, "ref": self.ref, "alt": self.alt})

        if len(r) > 0:
            func_class = r["cls"].squeeze()
        else:
            func_class = None

        return func_class

    def get_missplicing(self):

        if self.is_canonical:

            q = text(
                ''' SELECT * FROM vault WHERE chr = :chr AND pos = :pos AND ref = :ref AND alt LIKE :alt AND strand = :strand ''')

            r = pd.read_sql_query(q, connection,
                                  params={"chr": self.chr, "pos": self.pos, "ref": self.ref, "alt": f"%{self.alt}%",
                                          "strand": self.strand})

            if len(r) > 0:

                r = r.iloc[0]

                frame = r["frame"]

                if frame == "Start_Lost":
                    return ["STARTLOST", None, None, None, None, None, None]

                elif frame == "Stop_Lost":
                    return ["STOPLOST", None, None, None, None, None, None]

                elif frame == "UTR":
                    return ["UTR", None, None, None, None, None, None]

                elif frame == "F" or frame == "Inframe_w_stop_codon":

                    is_frameshifting = True

                elif frame == "I":
                    is_frameshifting = False

                trunc_start = r["trunc_start"]
                trunc_end = r["trunc_end"]
                csq_imp = r["csq_imp"]
                perc_lost = r["perc_lost"]
                NMD_exception = r["escapes_nmd"]

                if frame == "Inframe_w_stop_codon":
                    csq_imp = "STOPGAINED"

                if csq_imp == "ES":
                    exons_skipped = r["csq_pos"].split("-")
                    exons_skipped = [int(x) for x in exons_skipped]
                else:
                    exons_skipped = None

                if NMD_exception == 1:
                    undergoes_NMD = False
                else:
                    undergoes_NMD = self.check_NMD(trunc_start)

                return [csq_imp, is_frameshifting, undergoes_NMD, perc_lost, trunc_start, trunc_end, exons_skipped]

            else:
                return [None, None, None, None, None, None, None]

        else:

            return [None, None, None, None, None, None, None]

    def get_lowest_freq(self, exons):

        #exon_ranks = ",".join(map(str, exons))
        exons = tuple(exons)

        if self.strand == 1:
            q_exon = text(
                " SELECT * FROM exons WHERE enst_id = :tcpt_id AND exon_rank IN :exon_ranks ORDER BY highest_freq ASC LIMIT 1 ")
        else:
            q_exon = text(
                " SELECT * FROM exons WHERE enst_id = :tcpt_id AND exon_rank IN :exon_ranks ORDER BY highest_freq ASC LIMIT 1 ")

        exon = pd.read_sql_query(q_exon, connection, params={"tcpt_id": self.tcpt_id, "exon_ranks": exons})

        highest_freq = exon["highest_freq"].squeeze()

        if not highest_freq:
            highest_freq = 0

        return highest_freq


class ACMG:
    """
    Class for ACMG variant classification.
    Stores criteria evaluation, descriptions, flags, and final classification.
    """

    known_risk_variants = ["3-165830741-T-C", "7-117509089-C-T", "1-169549811-C-T", "X-154536002-C-T",
                            "9-36217448-C-T", "6-26090951-C-G", "6-26092913-G-A"]

    def __init__(self, variant: "Variant", predictor: Any) -> None:
        """
        Initialize ACMG classification object for a given variant.

        Args:
            variant (Variant): A Variant object containing effect and frequency data.
            predictor (Any): Predictor or model used to evaluate ACMG criteria.
        """
        self.variant: "Variant" = variant
        self.predictor: Any = predictor

        # ACMG criteria: PVS, PS, PM, PP, BA, BS, BP
        self.criteria: Dict[str, Optional[str]] = {k: None for k in [
            "PVS1", "PS1", "PS2", "PS3", "PS4",
            "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
            "PP1", "PP2", "PP3", "PP4", "PP5",
            "BA1", "BS1", "BS2", "BS3", "BS4",
            "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"
        ]}

        # Optional descriptive info for each criterion
        self.crit_desc: Dict[str, Optional[str]] = {k: None for k in self.criteria.keys()}

        # Any flags or notes about classification
        self.flags: List[str] = []

        # Final ACMG classification (Pathogenic, Likely Pathogenic, VUS, etc.)
        self.classification: Optional[str] = None

    def calc_ACMG_class(self, deactivated_criteria = None):

        criteria_met = []
        odds_path = 1
        post_prob = 0.1

        strength_key = {"A": float('inf'), "VS": 8, "VS-": 6, "S": 4, "S-": 3, "M": 2, "P": 1}
        odds_path_key = {"A": float('inf'), "VS": 350, "VS-": 350 ** 0.75, "S": 350 ** 0.5, "S-": 350 ** 0.375, "M": 350 ** 0.25, "P": 350 ** 0.125}

        total_point = 0

        for crit, assn in self.criteria.items():
            if (deactivated_criteria is None or crit not in deactivated_criteria) and (assn is not None) and (assn != 0):

                criteria_met.append(assn)

                crt_split = assn.split("_")
                impact = crt_split[0]
                strength = crt_split[1]


                if impact[0] == "B":
                    k = -1
                    odds_path = odds_path * (odds_path_key[strength]**-1)
                else:
                    k = 1
                    odds_path = odds_path * odds_path_key[strength]

                post_prob = (odds_path*0.1)/((odds_path-1)*0.1+1)

                point = strength_key[strength] * k
                total_point += point

        if total_point <= -7:
            self.classification = "BENIGN"
        elif total_point <= -1:
            self.classification = "LIKELY_BENIGN"
        elif total_point <= 1:
            self.classification = "VUS-LOW"
        elif total_point <= 3:
            self.classification = "VUS-MID"
        elif total_point <= 5:
            self.classification = "VUS-HIGH"
        elif total_point <= 9:
            self.classification = "LIKELY_PATHOGENIC"
        else:
            self.classification = "PATHOGENIC"

        self.variant.ACMG_criteria = criteria_met
        self.variant.ACMG_flags = sorted(set(self.flags))
        self.variant.ACMG_classification = self.classification
        self.variant.ACMG_score = total_point
        self.variant.post_prob = f"{round(post_prob*100,2)}%"

    def classify(self, PS4=None, deactivate_PM2=False, deactivate_PP5_BP6=False):

        deactivated_criteria = []

        if deactivate_PM2:
            deactivated_criteria.append("PM2")

        if deactivate_PP5_BP6:
            deactivated_criteria.append("PP5")
            deactivated_criteria.append("BP6")

        self.PS1()
        self.PS2()
        self.PS3_BS3()
        # self.PS4()
        self.PM1()
        self.PM2()
        self.PM3()
        self.PM4_BP3()
        self.PM5()
        # self.PM6()
        # self.PP1()
        self.PP2()
        self.PP3_BP4()
        # self.PP4()
        self.PP5_BP6()
        if self.variant.var_id not in ACMG.known_risk_variants:
            self.BA1()
            self.BS1()
            self.BS2()
        else:
            self.flags.append("COMMON_KNOWN_RISK_VARIANT")
        # self.BS4()
        self.BP1()
        # self.BP2()
        # self.BP5
        self._PVS1()

        self.__BP7()

        if PS4 == "PS4_S":
            self.criteria["PS4"] = "PS4_S"
        elif PS4 == "PS4_M":
            self.criteria["PS4"] = "PS4_M"
        else:
            self.criteria["PS4"] = None

        self.calc_ACMG_class(deactivated_criteria)

    # pop freq rules

    def PM2(self):

        # get
        bs1_cutoff = self.variant.BS1_cutoff

        if bs1_cutoff == -1:
            if self.variant.filt_af == 0:
                self.criteria["PM2"] = "PM2_P"
                self.flags.append("RARE_VARIANT")

            elif self.variant.MOI in ["BOTH", "BIALLELIC", "AR"] and self.variant.filt_af <= 0.01:
                self.criteria["PM2"] = "PM2_P"
                self.flags.append("RARE_VARIANT")

            elif (not self.variant.MOI or self.variant.MOI in ["MONOALLELIC", "AD", "X-LINKED",
                                                               "XL"]) and self.variant.filt_af <= 0.005:
                self.criteria["PM2"] = "PM2_P"
                self.flags.append("RARE_VARIANT")

            else:
                self.criteria["PM2"] = 0

        else:
            if self.variant.filt_af == 0:
                self.criteria["PM2"] = "PM2_P"
                self.flags.append("RARE_VARIANT")
            elif self.variant.filt_af <= bs1_cutoff / 10:
                self.criteria["PM2"] = "PM2_P"
                self.flags.append("RARE_VARIANT")
            else:
                self.criteria["PM2"] = 0

    def BA1(self):

        self.criteria["BA1"] = 0

        is_exception = False

        BA1_exceptions = {"ENSG00000244734": 0.05479,   #HBB
                          "ENSG00000010704": 0.14959}   #HFE

        BA1_cutoff = 0.05

        if self.variant.gene_id in BA1_exceptions.keys():
            BA1_cutoff = BA1_exceptions[self.variant.gene_id]
            is_exception = True

        if self.variant.filt_af > BA1_cutoff:
            self.criteria["BA1"] = "BA1_A"
            self.flags.append("POLYMORPHISM")

        if is_exception and self.criteria["BA1"] != "BA1_A":
                self.flags.append("BA1_EXCEPTION")

    def BS1(self):

        mocopav = self.variant.BS1_cutoff

        if mocopav == -1:
            if self.variant.filt_af > 0.015:
                self.criteria["BS1"] = "BS1_S"
                self.flags.append("TOO_COMMON_FOR_DISEASE")

            else:
                self.criteria["BS1"] = 0

        else:
            if self.variant.filt_af > mocopav:
                self.criteria["BS1"] = "BS1_S"
                self.flags.append("TOO_COMMON_FOR_DISEASE")
            else:
                self.criteria["BS1"] = 0

    def BS2(self):

        bs2_cutoff = self.variant.BS2_cutoff

        if bs2_cutoff != -1:
            if self.variant.hom_ac > bs2_cutoff:
                self.criteria["BS2"] = "BS2_P"
                self.flags.append("OBSERVED_HOMOZYGOUS")
            else:
                self.criteria["BS2"] = 0

        else:
            if self.variant.hom_ac > 2:
                self.criteria["BS2"] = "BS2_P"
                self.flags.append("OBSERVED_HOMOZYGOUS")

            else:
                self.criteria["BS2"] = 0

    # reported variant rules

    def PS1(self):

        if "missense" in self.variant.major_conseq:

            q_ps1 = text(''' SELECT * FROM clinvar 
                                WHERE enst_id LIKE :tcpt 
                                AND is_missense = 1
                                AND ABS(pos - :pos) <= 2 
                                AND aa_change LIKE :aa_change 
                                AND nt_change NOT LIKE :nt_change 
                                AND (cls = 'P' OR cls = 'LP') ''')

            r_ps1 = pd.read_sql_query(q_ps1, connection, params={"tcpt": f"%{self.variant.tcpt_id}%",
                                                                 "pos": self.variant.pos,
                                                                 "aa_change": f"%{self.variant.aa_change}%",
                                                                 "nt_change": f"%{self.variant.nt_change}%"})

            if len(r_ps1) > 1:
                self.criteria["PS1"] = "PS1_S"
                self.flags.append("SAME_AA_CHANGE_REPORTED_P/LP")
            elif len(r_ps1) > 0:
                star = r_ps1["star"].squeeze()
                if star == 0:
                    self.criteria["PS1"] = "PS1_P"
                elif star == 1:
                    self.flags.append("SAME_AA_CHANGE_REPORTED_P/LP")
                    self.criteria["PS1"] = "PS1_M"
                elif star == 2:
                    self.flags.append("SAME_AA_CHANGE_REPORTED_P/LP")
                    self.criteria["PS1"] = "PS1_S"
                elif star >= 3:
                    self.flags.append("SAME_AA_CHANGE_REPORTED_P/LP")
                    self.criteria["PS1"] = "PS1_VS"
            else:
                self.criteria["PS1"] = 0

    def PS2(self):
        if self.variant.is_denovo == 1:
            self.criteria["PS2"] = "PS2_S"
            self.flags.append("DE_NOVO")
        else:
            self.criteria["PS2"] = 0

    def PS3_BS3(self):

        rule_adjustments = {"ENSG00000177000": {"PS3": "S", "BS3": "P"}, #MTHFR
                       "ENSG00000157978": {"PS3": "S", "BS3": "S"}, #LDLRAP1
                       "ENSG00000012048": {"PS3": "S", "BS3": "S"}, #BRCA1
                       "ENSG00000075043": {"PS3": "S"}, #KCNQ2
                       "ENSG00000160200": {"PS3": "M", "BS3": "M"}, #CBS
                       "ENSG00000134086": {"PS3": "S"}, #VHL
                       "ENSG00000163914": {"PS3": "S"}, #RHO
                       "ENSG00000101981": {"PS3": "S"}, #F9
                       "ENSG00000215301": {"PS3": "S"}, #DDX3X
                       "ENSG00000171862": {"PS3": "M", "BS3": "M"}, #PTEN
                       "ENSG00000256269": {"PS3": "M", "BS3": "S"}, #HMBS
                       "ENSG00000198668": {"PS3": "S", "BS3": "M"}, #CALM1
                       "ENSG00000141510": {"PS3": "M", "BS3": "M"}, #TP53
                       "ENSG00000095002": {"PS3": "S", "BS3": "M"}, #MSH2
                       "ENSG00000106633": {"PS3": "M", "BS3": "P"}, #GCK
                       "ENSG00000203879": {"BS3": "S"}, #GDI1
                       "ENSG00000163930": {"PS3": "S", "BS3": "S"}, #BAP1
                       "ENSG00000036473": {"PS3": "S", "BS3": "P"}, #OTC
                       }

        func_class = self.variant.get_func_class()

        if func_class is not None:

            PS3_strength = "S"
            BS3_strength = "M"

            if self.variant.gene_id in rule_adjustments.keys():
                gene_adj = rule_adjustments[self.variant.gene_id]
                if "PS3" in gene_adj:
                    PS3_strength = gene_adj["PS3"]
                if "BS3" in gene_adj:
                    BS3_strength = gene_adj["BS3"]

            if func_class == 0:

                self.flags.append("IN_VITRO_NEUTRAL_IMPACT")

                self.criteria["BS3"] = f"BS3_{BS3_strength}"
                self.criteria["PS3"] = 0

            elif func_class == 2:

                self.flags.append("IN_VITRO_DAMAGING_IMPACT")

                self.criteria["BS3"] = 0
                self.criteria["PS3"] = f"PS3_{PS3_strength}"

        else:

            self.criteria["BS3"] = 0
            self.criteria["PS3"] = 0

    def PM3(self):

        if self.variant.is_intrans == 1:
            self.criteria["PM3"] = "PM3_M"
            self.flags.append("IN_TRANS")
        else:
            self.criteria["PM3"] = 0

    def PM5(self):

        if "missense" in self.variant.major_conseq:

            q_pm5 = text(''' SELECT * FROM clinvar WHERE enst_id LIKE :tcpt AND is_missense = 1
                                        AND aa_change NOT LIKE :aa_change 
                                        AND ABS(pos - :pos) <= 2 
                                        AND CONCAT(';',aa_pos,';') LIKE CONCAT('%;',:aa_pos,';%')
                                        AND (cls = 'P' OR cls = 'LP')  ''')

            r_pm5 = pd.read_sql_query(q_pm5, connection,
                                      params={"tcpt": f"%{self.variant.tcpt_id}%",
                                              "aa_change": f"%{self.variant.aa_change}%",
                                              "pos": self.variant.pos,
                                              "aa_pos": self.variant.aa_pos})
            if len(r_pm5) > 1:

                self.criteria["PM5"] = "PM5_S"
                self.flags.append("DIFFERENT_CHANGES_IN_SAME_CODON_REPORTED_P/LP")

            elif len(r_pm5) > 0:

                star = r_pm5["star"].squeeze()
                if star == 0:
                    self.criteria["PM5"] = 0
                elif star == 1:
                    self.flags.append("DIFFERENT_CHANGE_IN_SAME_CODON_REPORTED_P/LP")
                    self.criteria["PM5"] = "PM5_P"
                elif star >= 2:
                    self.flags.append("DIFFERENT_CHANGE_IN_SAME_CODON_REPORTED_P/LP")
                    self.criteria["PM5"] = "PM5_M"

            else:
                self.criteria["PM5"] = 0

    def PP5_BP6(self):

        if len(self.variant.clinvar) < 1:
            self.criteria["PP5"] = 0
            self.criteria["BP6"] = 0
        else:
            star = self.variant.clinvar["star"].squeeze()
            cls = self.variant.clinvar["cls"].squeeze()

            if star == 4:
                strength = "A"
            elif star == 3:
                strength = "VS"
            elif star == 2:
                strength = "S"
            elif star == 1:
                strength = "M"
            else:
                strength = "P"

            if cls in ["P", "LP"]:

                self.criteria["PP5"] = f"PP5_{strength}"
                self.flags.append(f"REPORTED_P/LP_{star}_STAR(S)")

            elif cls in ["B", "LB"]:

                self.criteria["BP6"] = f"BP6_{strength}"
                self.flags.append(f"REPORTED_B/LB_{star}_STAR(S)")

    # indel & stop-loss rules

    def PM4_BP3(self):


        if self.variant.major_conseq == "stop_lost":
            if self.variant.alt_stop != -1:
                self.criteria["PM4"] = "PM4_M"
                self.flags.append(f"PROTEIN_LENGTH_CHANGED")
            else:
                self.criteria["PM4"] = 0

        elif self.variant.major_conseq in ["inframe_insertion", "inframe_deletion"]:
            #print("inframe indel")
            if self.variant.check_repeat()["within_repeat"] is True:  # checks if within repeat region
                #print("within repeat")
                if self.variant.check_repeat()[
                    "contains_pathogenic"] is False:  # if yes, checks if region absent from pat vars

                    self.flags.append(f"IN_REDUNDANT_REPEAT")
                    self.criteria["BP3"] = "BP3_P"

                else:
                    self.criteria["PM4"] = "PM4_M"
                    self.flags.append("CRITICAL_REPEAT")
            else:
                self.criteria["PM4"] = "PM4_M"
                self.flags.append("PROTEIN_LENGTH_CHANGES_IN_NON-REPEAT_REGION")	
        else:
            self.criteria["PM4"] = 0

    # in-silico prediction

    def PP3_BP4(self):
        if self.variant.is_splicing is True and isinstance(self.variant.spliceAI, float):

            if self.variant.spliceAI >= 0.2:
                self.flags.append("PATHOGENIC_IN_SILICO_PREDICTION")
                self.criteria["PP3"] = "PP3_M"
                self.criteria["BP4"] = 0
            elif self.variant.spliceAI <= 0.1:
                self.flags.append("BENIGN_IN_SILICO_PREDICTION")
                self.criteria["PP3"] = 0
                self.criteria["BP4"] = "BP4_M"
        else:

            if self.predictor == "revel":

                if not self.variant.REVEL:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = 0
                elif self.variant.REVEL >= 0.932:
                    self.criteria["PP3"] = "PP3_S"
                    self.criteria["BP4"] = 0
                    self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                elif self.variant.REVEL >= 0.773:
                    self.criteria["PP3"] = "PP3_M"
                    self.criteria["BP4"] = 0
                    self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                elif self.variant.REVEL >= 0.644:
                    self.criteria["PP3"] = "PP3_P"
                    self.criteria["BP4"] = 0
                    self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                elif self.variant.REVEL > 0.290:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = 0
                elif self.variant.REVEL > 0.183:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = "BP4_P"
                    self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")
                elif self.variant.REVEL > 0.016:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = "BP4_M"
                    self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")
                elif self.variant.REVEL >= 0:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = "BP4_S"
                    self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")

            if self.predictor == "bayesdel":

                if not self.variant.BayesDel:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = 0
                elif self.variant.BayesDel >= 0.5:
                    self.criteria["PP3"] = "PP3_S"
                    self.criteria["BP4"] = 0
                    self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                elif self.variant.BayesDel >= 0.27:
                    self.criteria["PP3"] = "PP3_M"
                    self.criteria["BP4"] = 0
                    self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                elif self.variant.BayesDel >= 0.13:
                    self.criteria["PP3"] = "PP3_P"
                    self.criteria["BP4"] = 0
                    self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                elif self.variant.BayesDel > -0.18:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = 0
                elif self.variant.BayesDel > -0.36:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = "BP4_P"
                    self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")
                elif self.variant.BayesDel > -1.5:
                    self.criteria["PP3"] = 0
                    self.criteria["BP4"] = "BP4_M"
                    self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")

    # hotspot rule

    def PM1(self):

        if "missense" in self.variant.major_conseq:

            if self.variant.domain:

                clv_stats = self.variant.clinvar_within(method="domain")
                pat_vars, ben_vars = clv_stats

                if pat_vars / (pat_vars + ben_vars + 0.5) > 0.5:
                    self.criteria["PM1"] = "PM1_M"
                    self.flags.append(f"CRITICAL_DOMAIN")
                else:
                    self.criteria["PM1"] = 0

            elif self.variant.check_regional_constraint(): # and self.variant.get_pext() and self.variant.get_pext() > 0.5:

                self.criteria["PM1"] = "PM1_M"
                self.flags.append(f"CONSTRAINED_REGION")

            else:
                clv_stats = self.variant.clinvar_within(method="exon")
                pat_vars, ben_vars = clv_stats
                if pat_vars / (pat_vars + ben_vars + 0.5) > 0.5:  # ben_vars == 0 and pat_vars >= 2:
                    self.criteria["PM1"] = "PM1_M"
                    self.flags.append(f"EXONIC_HOTSPOT")
                else:
                    self.criteria["PM1"] = 0

        else:
            self.criteria["PM1"] = 0

    # gene missense sensitivity rules

    def PP2(self):

        if "missense" in self.variant.major_conseq:

            if self.variant.mis_z > 3.09:
                self.criteria["PP2"] = "PP2_P"
                self.flags.append(f"MISSENSE_INTOLERANT_GENE")

            else:
                clv = self.variant.clinvar_within(method="gene")

                pat_vars, ben_vars = clv
                if pat_vars / (ben_vars + pat_vars + 0.5) > 0.5:
                    self.criteria["PP2"] = "PP2_P"
                    self.flags.append(f"ENRICHED_FOR_P/LP_MISSENSE")

                else:
                    self.criteria["PP2"] = 0
        else:
            self.criteria["PP2"] = 0

    def BP1(self):

        if "missense" in self.variant.major_conseq:
            #print(self.variant.lof_perc)
            if self.variant.lof_perc > 0.85:
                self.criteria["BP1"] = "BP1_P"
                self.flags.append(f"MISSENSE_IN_LOF_GENE")

            else:
                self.criteria["BP1"] = 0

        else:
            self.criteria["BP1"] = 0

    # loss-of-function rule
    def _PVS1(self, ignore_lof_check=True):

        is_null = False
        for conseq in self.variant.all_conseqs:
            if "stop_gained" == conseq:
                is_null = True
                break
            elif "stop_lost" == conseq:
                is_null = True
                break
            elif "frameshift" in conseq:
                is_null = True
                break
            elif conseq in ["start_lost", "initiator_codon_variant"]:
                is_null = True
                break

        # for some reason, all_conseq and major_conseq from ensembl dont match sometimes
        if self.variant.aa_change:
            if self.variant.major_conseq in ["stop_gained", "frameshift_variant", "start_lost", "stop_lost"
                                                                                                "initiator_codon_variant"]:
                is_null = True

        if is_null:
            if ignore_lof_check is True or self.variant.is_in_LOF_gene() is True:

                if ignore_lof_check is False:
                    self.flags.append("LOF_INTOLERANT_GENE")

                if self.variant.major_conseq in ["stop_gained",
                                                 "frameshift_variant"] and "stop_retained_variant" not in self.variant.all_conseqs:

                    if self.variant.major_conseq == "frameshift_variant":
                        if self.variant.check_homopolymer():
                            self.flags.append("FRAMESHIFT_IN_HOMOPOLYMER")

                        if self.variant.check_repeat()["within_repeat"]:
                            self.flags.append("FRAMESHIFT_IN_REPEAT")

                    #pext = self.variant.get_pext()

                    if self.variant.undergoes_NMD:

                        self.flags.append("UNDERGOES_NMD")

                        #if pext is None or (pext and pext > 0.5):

                        self.criteria["PVS1"] = "PVS1_VS"
                        #self.flags.append("HIGH_PEXT")

                        #else:

                         #   self.criteria["PVS1"] = 0
                          #  self.flags.append("LOW_PEXT")

                    else:

                        self.flags.append("ESCAPES_NMD")

                        clv_stats = self.variant.clinvar_within()

                        pat_vars, ben_vars = clv_stats
                        tot_vars = pat_vars + ben_vars + 0.5

                        if pat_vars / tot_vars > 0.5:
                            self.criteria["PVS1"] = "PVS1_S"
                            self.flags.append("REMOVES_EXON_ENRICHED_FOR_P/LP_VARIANTS")

                        elif self.variant.check_regional_constraint():

                            self.criteria["PVS1"] = "PVS1_S"
                            self.flags.append("CONSTRAINED_REGION")

                        elif self.variant.check_domains(extent="tcpt") and self.variant.is_domain_lost():

                            self.criteria["PVS1"] = "PVS1_S"
                            self.flags.append("ALTERS_FUNCTIONAL_DOMAIN")

                        else:

                            #if not pext or (pext and pext > 0.5):

                            #self.flags.append("HIGH_PEXT")

                            if self.variant.exon_start:
                                highest_af = self.variant.highest_freq
                            else:
                                highest_af = 0

                            if highest_af >= 0.01:
                                self.criteria["PVS1"] = 0
                                self.flags.append(f"EXON_HARBORS_COMMON_LOF")

                            else:

                                # in case, the frameshift variant lies at exon-intron junction
                                if not self.variant.aa_pos:
                                    if "stop_lost" in self.variant.all_conseqs:
                                        closest_aa = self.variant.prot_length
                                    else:
                                        closest_aa = math.ceil(
                                            int(re.findall(r"(\d+)(?=[+-_])", self.variant.nt_change)[0]) / 3)
                                else:
                                    closest_aa = self.variant.aa_pos

                                perc_lost = 1 - (closest_aa / self.variant.prot_length)

                                

                                if perc_lost >= 0.1:
                                    self.criteria["PVS1"] = "PVS1_S"
                                    self.flags.append("AT_LEAST_10%_OF_PROTEIN_LOST")

                                else:
                                    self.criteria["PVS1"] = "PVS1_M"
                                    self.flags.append("LESS_THAN_10%_OF_PROTEIN_LOST")

                elif self.variant.major_conseq in ["stop_lost"]:

                    if self.variant.alt_stop == -1:
                        self.criteria["PVS1"] = "PVS1_VS"
                        self.flags.append("NO_ALTERNATIVE_STOP")
                    elif self.variant.alt_stop == -2:
                        self.criteria["PVS1"] = 0
                        self.flags.append("UNDEFINED_TCPT_SEQ")
                    else:
                        self.criteria["PVS1"] = 0

                elif self.variant.major_conseq in ["start_lost"]:

                    if self.variant.alt_start > 0:
                        self.flags.append("RESCUE_START_CODON")

                        perc_alt_start = self.variant.alt_start / self.variant.prot_length

                        if perc_alt_start <= 0.25:

                            self.flags.append("LESS_THAN_25%_OF_PROTEIN_LOST")

                            if self.variant.strand == 1:
                                pos_alt_start = self.variant.tln_start + self.variant.alt_start * 3 - 2
                            else:
                                pos_alt_start = self.variant.tln_start - self.variant.alt_start * 3 + 2

                            pat_vars = self.variant.clinvar_within(method="alternative_start")[0]
                            domains = self.variant.check_domains(extent="tcpt")

                            domain_lost = False

                            if domains:

                                for domain, domain_range in domains.items():

                                    domain_start = domain_range[0]
                                    domain_end = domain_range[1]

                                    if self.variant.strand == 1:

                                        if (self.variant.tln_start <= domain_start <= pos_alt_start) or (
                                                self.variant.tln_start <= domain_end <= pos_alt_start):
                                            domain_lost = True

                                    else:

                                        if (self.variant.tln_start >= domain_start >= pos_alt_start) or (
                                                self.variant.tln_start >= domain_end >= pos_alt_start):
                                            domain_lost = True

                            if pat_vars >= 1:
                                self.flags.append("AT_LEAST_ONE_P/LP_REPORTED_IN_LOST_REGION")

                            if pat_vars >= 1 or domain_lost:
                                self.criteria["PVS1"] = "PVS1_M"

                            else:
                                self.criteria["PVS1"] = "PVS1_P"

                        else:
                            self.flags.append("MORE_THAN_25%_OF_PROTEIN_LOST")
                            self.criteria["PVS1"] = "PVS1_M"
                    else:
                        self.criteria["PVS1"] = "PVS1_M"
                        self.flags.append("NO_RESCUE_START_CODON_IN_TRANSCRIPT")

                else:
                    self.criteria["PVS1"] = 0
            else:
                self.flags.append("LOF_TOLERANT_GENE")
                self.criteria["PVS1"] = 0
        else:
            self.criteria["PVS1"] = 0

    def __BP7(self):

        if self.variant.is_splicing is True:

            annotated_splicing = self.variant.get_missplicing()

            if annotated_splicing[0] is not None:
                spl_csq, is_frameshifting, undergoes_NMD, perc_lost, trunc_start, trunc_end, exons_skipped = annotated_splicing
		
                #print(f"percentage lost: {round(perc_lost,2)*100}%")
		
                if is_frameshifting is True:

                    if spl_csq == "STOPGAINED":
                        self.flags.append("RETAINED_INTRON_W_STOP_CODON")
                        self.flags.append("INFRAME")

                        self.variant.all_conseqs.append("stop_gained")
                    else:
                        self.flags.append("FRAMESHIFT")

                    if undergoes_NMD:
                        self.flags.append("UNDERGOES_NMD")
                        #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                        self.criteria["PVS1"] = "PVS1_VS-"
                    else:
                        self.flags.append("ESCAPES_NMD")

                        clv_stats = self.variant.clinvar_within(method="range", query_start=trunc_start,
                                                                query_end=trunc_end)

                        pat_vars, ben_vars = clv_stats
                        tot_vars = pat_vars + ben_vars + 0.5

                        if pat_vars / tot_vars > 0.5:
                            self.criteria["PVS1"] = "PVS1_S-"
                            #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                            self.flags.append("EXON_ENRICHED_FOR_P/LP_VARIANTS")

                        elif self.variant.check_regional_constraint(trunc_start, trunc_end):

                            self.criteria["PVS1"] = "PVS1_S-"
                            #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                            self.flags.append("CONSTRAINED_REGION")

                        elif self.variant.check_domains(extent="range", query_start=trunc_start, query_end=trunc_end):
                            self.criteria["PVS1"] = "PVS1_S-"
                            #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                            self.flags.append("ALTERS_FUNCTIONAL_DOMAIN")

                        else:

                            if self.variant.exon_start:
                                highest_af = self.variant.highest_freq
                            else:
                                highest_af = 0

                            if highest_af >= 0.01:
                                self.criteria["PVS1"] = 0
                                self.flags.append("EXON_HARBORS_COMMON_LOF")

                            else:

                                if perc_lost >= 0.1:
                                    self.criteria["PVS1"] = "PVS1_S-"
                                    #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                                    self.flags.append("AT_LEAST_10%_OF_PROTEIN_LOST")

                                else:
                                    self.criteria["PVS1"] = "PVS1_M"
                                    self.flags.append("LESS_THAN_10%_OF_PROTEIN_LOST")

                elif spl_csq == "STARTLOST":

                    self.flags.append("MIS-SPLICING_ABOLISHING_START_CODON")
                    self.criteria["PVS1"] = "PVS1_P"

                elif spl_csq == "STOPLOST":

                    self.flags.append("MIS-SPLICING_ABOLISHING_STOP_CODON")
                    self.criteria["PVS1"] = "PM4_M"

                elif spl_csq == "UTR":
                    self.flags.append("SPLICE_AFFECTING_UTR")

                elif is_frameshifting is False:

                    self.flags.append("INFRAME")

                    clv_stats = self.variant.clinvar_within(method="range", query_start=trunc_start,
                                                            query_end=trunc_end)

                    pat_vars, ben_vars = clv_stats
                    tot_vars = pat_vars + ben_vars + 0.5

                    if self.variant.check_domains(extent="range", query_start=trunc_start, query_end=trunc_end):
                        self.criteria["PVS1"] = "PVS1_S-"
                        #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                        self.flags.append("ALTERS_FUNCTIONAL_DOMAIN")

                    elif pat_vars / tot_vars > 0.5:

                        self.criteria["PVS1"] = "PVS1_S-"
                        #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                        self.flags.append("EXON_ENRICHED_FOR_P/LP_VARIANTS")

                    elif self.variant.check_regional_constraint(trunc_start, trunc_end):

                        self.criteria["PVS1"] = "PVS1_S-"
                        #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                        self.flags.append("CONSTRAINED_REGION")

                    else:

                        if spl_csq == "ES":
                            highest_af = self.variant.get_lowest_freq(exons_skipped)
                        else:
                            highest_af = self.variant.highest_freq  # self.variant.lof_within_exon()

                        if highest_af >= 0.01:
                            self.flags.append("LOF_VARIANTS_COMMON_IN_EXON")
                            self.criteria["PVS1"] = 0

                        else:

                            if perc_lost >= 0.1:
                                self.criteria["PVS1"] = "PVS1_S-"
                                #self.flags.append("COULD_BE_RESCUED_BY_ALTERNATIVE_SPLICING")
                                self.flags.append("AT_LEAST_10%_OF_PROTEIN_LOST")
                            else:
                                self.flags.append("LESS_THAN_10%_OF_PROTEIN_LOST")
                                self.criteria["PVS1"] = "PVS1_M"

                            ###from BP7
                            clvspl = self.variant.clinvar_within(method="splice_site")
                            if clvspl:
                                if clvspl[0] > 0 or clvspl[4] > 0:
                                    self.criteria["PS1"] = "PS1_S"
                                elif clvspl[2] > 3:
                                    self.criteria["PS1"] = "PS1_M"
                                elif clvspl[3] > 0:
                                    self.criteria["PS1"] = "PS1_P"

                clvspl = self.variant.clinvar_within(method="splice_site")

                if clvspl:

                    if clvspl[0] > 0 or clvspl[4] > 0:

                        self.criteria["PS1"] = "PS1_S"

                    elif clvspl[2] > 3:

                        self.criteria["PS1"] = "PS1_M"

                    elif clvspl[3] > 0:

                        self.criteria["PS1"] = "PS1_P"

            else:

                if self.variant.is_at_canonical_splice_site():

                    clvspl = self.variant.clinvar_within(method="splice_site")

                    if clvspl:

                        if clvspl[0] > 0 or clvspl[4] > 0:

                            self.criteria["PS1"] = "PS1_S"

                        elif clvspl[2] > 3:

                            self.criteria["PS1"] = "PS1_M"

                        elif clvspl[3] > 0:

                            self.criteria["PS1"] = "PS1_P"

                else:

                    if self.variant.spliceAI:
                        if self.variant.spliceAI >= 0.2:
                            self.flags.append("PATHOGENIC_IN_SILICO_PREDICTION")
                            self.criteria["PP3"] = "PP3_M"
                            clvs = self.variant.clinvar_within(method="splice_site")

                            if not clvs:
                                self.criteria["PS1"] = 0
                            elif clvs[0] > 0:
                                self.criteria["PS1"] = "PS1_S"
                            elif clvs[1] > 0:
                                self.criteria["PS1"] = "PS1_M"
                            elif clvs[2] > 0:
                                self.criteria["PS1"] = "PS1_M"
                            elif clvs[3] > 0:
                                self.criteria["PS1"] = "PS1_P"

                        elif self.variant.spliceAI <= 0.1:

                            if any(("synonymous" in effect) or ("intron" in effect) for effect in
                                   self.variant.all_conseqs):
                                self.criteria["BP4"] = "BP4_M"
                                self.criteria["BP7"] = "BP7_P"
                            elif all("missense" not in effect for effect in self.variant.all_conseqs):
                                self.criteria["BP4"] = "BP4_M"
                                self.criteria["PP3"] = 0
                            else:
                                # from PP3
                                if self.predictor == "revel":

                                    if not self.variant.REVEL:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = 0
                                    elif self.variant.REVEL >= 0.932:
                                        self.criteria["PP3"] = "PP3_S"
                                        self.criteria["BP4"] = 0
                                        self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                                    elif self.variant.REVEL >= 0.773:
                                        self.criteria["PP3"] = "PP3_M"
                                        self.criteria["BP4"] = 0
                                        self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                                    elif self.variant.REVEL >= 0.644:
                                        self.criteria["PP3"] = "PP3_P"
                                        self.criteria["BP4"] = 0
                                        self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                                    elif self.variant.REVEL > 0.290:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = 0
                                    elif self.variant.REVEL > 0.183:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = "BP4_P"
                                        self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")
                                    elif self.variant.REVEL > 0.016:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = "BP4_M"
                                        self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")
                                    elif self.variant.REVEL >= 0:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = "BP4_S"
                                        self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")

                                if self.predictor == "bayesdel":

                                    if not self.variant.BayesDel:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = 0
                                    elif self.variant.BayesDel >= 0.5:
                                        self.criteria["PP3"] = "PP3_S"
                                        self.criteria["BP4"] = 0
                                        self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                                    elif self.variant.BayesDel >= 0.27:
                                        self.criteria["PP3"] = "PP3_M"
                                        self.criteria["BP4"] = 0
                                        self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                                    elif self.variant.BayesDel >= 0.13:
                                        self.criteria["PP3"] = "PP3_P"
                                        self.criteria["BP4"] = 0
                                        self.flags.append(f"PATHOGENIC_IN_SILICO_PREDICTION")
                                    elif self.variant.BayesDel > -0.18:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = 0
                                    elif self.variant.BayesDel > -0.36:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = "BP4_P"
                                        self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")
                                    elif self.variant.BayesDel > -1.5:
                                        self.criteria["PP3"] = 0
                                        self.criteria["BP4"] = "BP4_M"
                                        self.flags.append(f"BENIGN_IN_SILICO_PREDICTION")

                        else:
                            pass  # consider protein-impact

        else:
            for conseq in self.variant.all_conseqs:
                if "synon" in conseq or "stream" in conseq or "prime" in conseq or "intron" in self.variant.major_conseq:
                    if self.variant.phyloP100way and self.variant.phyloP100way < 0.1:
                        self.criteria["BP7"] = "BP7_P"


class AAVC:
    """
    Automated ACMG-based Variant Classifier (AAVC) main class.
    Manages variant input, processing, and output generation.
    """

    def __init__(self,
                 variants: Optional[List[str]] = None,
                 mode: str = "default",
                 clock: bool = False,
                 debug: bool = False,
                 cache: bool = False,
                 prioritization: str = "canonical",
                 predictor: str = "bayesdel",
                 deactivate_PM2: bool = False,
                 deactivate_PP5_BP6: bool = False) -> None:
        """
        Initialize AAVC object.

        Args:
            variants (Optional[List[str]], optional): List of variant strings. Defaults to None.
            mode (str, optional): Processing mode. Defaults to "default".
            clock (bool, optional): Enable runtime tracking. Defaults to False.
            debug (bool, optional): Enable debug mode. Defaults to False.
            cache (bool, optional): Enable caching of intermediate results. Defaults to False.
            prioritization (str, optional): Variant prioritization strategy. Defaults to "canonical".
            predictor (str, optional): Predictor type ("bayesdel" or "revel"). Defaults to "bayesdel".
            deactivate_PM2 (bool, optional): Option to deactivate PM2 criterion. Defaults to False.
            deactivate_PP5_BP6 (bool, optional): Option to deactivate PP5/BP6 criteria. Defaults to False.
        """
        self.variants: List[str] = variants if variants else []
        self.output: Dict[str, Any] = {}
        self.var_objects: List[Any] = []
        self.mode: str = mode
        self.debug: bool = debug
        self.cache: bool = cache
        self.clock: bool = clock or debug  # runtime clock enabled if debug
        self.prioritization: str = prioritization
        self.predictor: str = predictor
        self.deactivate_PM2: bool = deactivate_PM2
        self.deactivate_PP5_BP6: bool = deactivate_PP5_BP6
        self.raised_exception: bool = False
        self.exception_statement: str = ""

    def process(self, variant):

        try:
            variant.classify(self.mode, self.deactivate_PM2, self.deactivate_PP5_BP6)

            print(f"\033[1mvariant: {variant.var_id} ({variant.gene}:{variant.nt_change}) ({variant.post_prob}) \033[1;97m<<{variant.ACMG_classification}>>\033[0m")
            print(variant.ACMG_criteria)
            print(variant.ACMG_flags)
            print(variant.eff_data)
            print(variant.freq_data)
            print("\n")

            self.output[variant.var_id] = {"Gene": variant.gene, "GeneID": variant.gene_id, "TcptID": variant.tcpt_id,
                                           "nt_change": variant.nt_change, "aa_change": variant.aa_change,
                                           "effect": "&".join(variant.all_conseqs),
                                           "REVEL": variant.REVEL, "BayesDel": variant.BayesDel,
                                           "SpliceAI": variant.spliceAI,
                                           "PhyloP100way": variant.phyloP100way, "gnomAD4_popmax": variant.filt_af,
                                           "gnomAD4_nhomalt": variant.hom_ac,
                                           "ACMG_Class": variant.ACMG_classification, "ACMG_Score": variant.ACMG_score,
                                           "ACMG_Criteria": variant.ACMG_criteria, "ACMG_Flags": variant.ACMG_flags}
            self.var_objects.append(variant)

        except Exception as e:
            self.raised_exception = True
            self.exception_statement = "Could not classify!"
            print(f"Could not classify {variant.var_id}: {e}")
            self.output[variant.var_id] = {"Gene": None, "GeneID": None, "TcptID": None,
                                           "nt_change": None, "aa_change": None,
                                           "effect": None,
                                           "REVEL": None, "BayesDel": None, "SpliceAI": None,
                                           "PhyloP100way": None, "gnomAD4_popmax": None,
                                           "gnomAD4_nhomalt": None,
                                           "ACMG_Class": None, "ACMG_Score": None,
                                           "ACMG_Criteria": None, "ACMG_Flags": None}
            self.var_objects.append(Variant)
            if self.debug:
                traceback.print_exc()

    def run(self, batch_size=25):

        try:

            clock_start = time.time()

            var_count = len(self.variants)

            # max 200 vars -- rate-limited at 55,000 per hour
            if var_count > batch_size:

                groups = []

                for i in range(0, var_count, batch_size):
                    group = self.variants[i:i + batch_size]
                    groups.append(group)

                for group in groups:
                    batch = Query(group, debug=self.debug, predictor=self.predictor)

                    for var_id in batch.eff_data.keys():
                        try:
                            eff = batch.eff_data[var_id]
                            freq = batch.freq_data[var_id]
                            # moi = batch.MOI_data[var_id]

                            variant = Variant(var_id, eff, freq, predictor=self.predictor)
                            self.process(variant)

                        except Exception as e:
                            self.raised_exception = True
                            self.exception_statement = "Query failed!"
                            print(f"Query failed for variant {var_id}:", str(e))
                            if self.debug:
                                traceback.print_exc()

            else:
                if isinstance(self.variants, list):
                    batch = Query(self.variants, debug=self.debug, cache=self.cache, predictor=self.predictor)

                else:
                    batch = Query([self.variants], debug=self.debug, cache=self.cache, predictor=self.predictor)

                for var_id in batch.eff_data.keys():
                    eff = batch.eff_data[var_id]
                    freq = batch.freq_data[var_id]
                    # moi = batch.MOI_data[var_id]

                    variant = Variant(var_id, eff, freq, predictor=self.predictor)

                    self.process(variant)

            if not self.debug:
                self.export()

            if self.clock:
                clock_end = time.time()
                print(f"elapsed: {round((clock_end - clock_start) / var_count, 2)} sec per variant")
                print(f"interpretation completed: {time.ctime()}")

        except Exception as e:
            self.raised_exception = True
            self.exception_statement = "Unknown error!"
            print("An unexpected error occurred:", str(e))
            if self.debug:
                traceback.print_exc()
            else:
                self.export()

    def vcf_run(self, vcf_dir, cohort_dir=None, mode="OB", svcf_mode=False):

        clock_start = time.time()

        if not svcf_mode:
            print("Parsing VCF...")
            # command = ["Rscript", os.path.abspath("sumvar2.R"), vcf_dir, cohort_dir, mode]
            command = ["Rscript", os.path.abspath("vcf_parser.R"), "--vcf_dir", vcf_dir, "--keep_info"]

            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

        o_name = vcf_dir.split('.')[0]

        if svcf_mode is True or process.returncode == 0:

            svcf = pd.read_csv(f"{o_name}.svcf", sep="\t", low_memory=False)
            svcf = svcf[svcf["nt_change"] != "."]
            svcf['id'] = svcf.iloc[:, :4].apply(lambda row: '-'.join(map(str, row)), axis=1)

            self.variants = svcf["id"].tolist()

            var_count = len(self.variants)

            if var_count == 0:
                var_count = 1

            for var in self.variants:

                var_row = svcf[svcf["id"] == var]

                gene = var_row["gene"].squeeze()
                gene_id = var_row["ensg_id"].squeeze()
                tcpt_id = var_row["enst_id"].squeeze()

                effect = var_row["effect"].squeeze()
                if "&" in effect:
                    major_conseq = effect.split("&")[0]
                    conseqs = effect.split("&")
                else:
                    major_conseq = effect
                    conseqs = [effect]

                nt_change = f"c.{Query.convert_nt_notation(var_row['nt_change'].squeeze())}"

                raw_aa_change = var_row["aa_change"].squeeze()

                if raw_aa_change == "." or raw_aa_change == "NA":
                    aa_change = None
                    aa_pos = None
                else:
                    aa_change = f"p.{Query.convert_aa_notation(raw_aa_change)}"
                    pos_match = re.match(r"p.[A-Za-z]([0-9]+).*", aa_change)
                    if pos_match and "deletion" not in effect and "insertion" not in effect:
                        aa_pos = int(pos_match.group(1))
                    else:
                        aa_pos = None

                mss_pred = var_row["predictor"].squeeze()
                mss_pred = None if mss_pred == "." or mss_pred == "NA" else float(mss_pred)

                revel, bayesdel = None, None
                if self.predictor == "revel":
                    revel = mss_pred
                elif self.predictor == "bayesdel":
                    bayesdel = mss_pred
                else:
                    print("Invalid missense predictor!")

                spliceai = var_row["spliceai"].squeeze()
                spliceai = None if spliceai == "." or spliceai == "NA" else float(spliceai)

                phylop100way = var_row["phylop"].squeeze()
                phylop100way = None if phylop100way == "." or phylop100way == "NA" else float(phylop100way)

                popmax = var_row["grpmax"].squeeze()
                popmax = 0 if popmax == "." or popmax == "NA" else float(popmax)

                PS4 = None

                total_hom_ac = var_row["nhomalt"].squeeze()
                total_hom_ac = 0 if total_hom_ac == "." or pd.isnull(total_hom_ac) or total_hom_ac == "NA" else int(
                    total_hom_ac)

                try:

                    q = text(" SELECT * FROM exons WHERE enst_id = :tcpt_id ")

                    gene_match = pd.read_sql_query(q, connection, params={"tcpt_id": tcpt_id})

                    if len(gene_match) > 0:
                        strand = gene_match["strand"][0].squeeze()
                        is_canonical = gene_match["is_canonical"][0].squeeze()
                        '''
                        else:
    
                            gene_query = text(" SELECT * FROM exons WHERE gene = :gene ")
                            gene_match = pd.read_sql_query(gene_query, connection, params={"gene": gene})
    
                            if len(gene_match) > 0:
                                strand = gene_match["strand"][0].squeeze()
                                print(strand)
                                is_canonical = 0
                        '''
                    else:
                        self.raised_exception = True
                        self.exception_statement = "Transcript could not be matched!"
                        raise Exception("Transcript could not be matched!")

                    var_data = {"gene": gene, "gene_id": gene_id, "tcpt_id": tcpt_id, "strand": strand,
                                "nt_change": nt_change, "aa_change": aa_change, "aa_pos": aa_pos,
                                "major_conseq": major_conseq, "bayesdel": bayesdel,
                                "all_conseqs": conseqs, "revel": revel, "spliceAI": spliceai,
                                "phyloP100way": phylop100way,
                                "is_canonical": is_canonical}

                    vua = Variant(var, var_data, {"popmax": popmax, "total_hom_ac": total_hom_ac},
                                  predictor=self.predictor)
                    vua.classify(PS4)
                    print("\n")
                    print(f"\033[1mvariant: {var} ({vua.gene}:{vua.nt_change}) \033[1;97m<<{vua.ACMG_classification}>>\033[0m")
                    print(vua.ACMG_criteria)
                    print(vua.ACMG_flags)
                    svcf.loc[svcf["id"] == var, "ACMG_class"] = vua.ACMG_classification
                    svcf.loc[svcf["id"] == var, "ACMG_score"] = vua.ACMG_score
                    svcf.loc[svcf["id"] == var, "ACMG_criteria"] = str(vua.ACMG_criteria)
                    svcf.loc[svcf["id"] == var, "ACMG_flags"] = str(vua.ACMG_flags)

                except Exception as e:
                    self.raised_exception = True
                    self.exception_statement = "Unknown error!"
                    print(f"Variant {var} could not be classified:", str(e))
                    traceback.print_exc()
                    svcf.loc[svcf["id"] == var, "ACMG_class"] = None
                    svcf.loc[svcf["id"] == var, "ACMG_score"] = None
                    svcf.loc[svcf["id"] == var, "ACMG_criteria"] = None
                    svcf.loc[svcf["id"] == var, "ACMG_flags"] = None

            filename = f"AAVC_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

            svcf.to_csv(f"{filename}", sep="\t", index=False)

            if self.clock:
                clock_end = time.time()
                print(f"elapsed: {round((clock_end - clock_start) / var_count, 2)} sec per variant")
                print(f"interpretation completed: {time.ctime()}")

        else:
            print(f"Error running R script: {stderr.decode('utf-8')}")

    def export(self):

        filename = f"AAVC_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        out_df = pd.DataFrame.from_dict(self.output, orient='index')

        out_df.to_csv(filename, sep="\t", index=True)

    def count_variants(self):

        result = subprocess.run(['wc', '-l', self.variants], capture_output=True, text=True)
        line_count = int(result.stdout.split()[0])
        return line_count

# uncomment below to run AAVC in Python interpreter

'''
# example variant list
variant_list = pd.read_csv("var.txt", sep="\t", header=None)[0].tolist()

# or run a single variant
variant_list = ["3-10142140-A-G"]

# create AAVC instance
run = AAVC(variant_list, clock=True, debug=True, predictor="bayesdel")

# run the classifier
run.run()
'''

if __name__ == "__main__":
    # check PostgreSQL status
    check_db = subprocess.run(
        ['systemctl', 'status', 'postgresql@16-main'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if 'active (running)' not in check_db.stdout:
        subprocess.run("sudo systemctl start postgresql", shell=True, check=True)

    # parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Ex: python aavc.py input.txt --vcf_mode --debug'
    )

    parser.add_argument('input', type=str, help='variant list (input.txt) or variant ID (12-106992962-T-G)')
    parser.add_argument('--debug', action='store_true', help='run the interpreter in debug mode')
    parser.add_argument('--vcf_mode', action='store_true', help='process a pre-calculated VCF file')
    parser.add_argument('--svcf_mode', action='store_true', help='process a pre-calculated sVCF file')
    parser.add_argument('--keep_info', action='store_true', help='keep the INFO column if the input is VCF')
    parser.add_argument('--predictor', choices=['bayesdel', 'revel'], help='missense prediction tool to use')
    parser.add_argument('--activate_PM2', action='store_true', help='activate rare variant criteria')
    parser.add_argument('--activate_PP5_BP6', action='store_true', help='activate reputable source criteria')

    args = parser.parse_args()

    PM2 = args.activate_PM2
    PP5_BP6 = args.activate_PP5_BP6

    run = None

    if args.vcf_mode:
        run = AAVC(clock=True, deactivate_PM2=not PM2, deactivate_PP5_BP6=not PP5_BP6, predictor=args.predictor)
        run.vcf_run(args.input)

    elif args.svcf_mode:
        run = AAVC(clock=True, deactivate_PM2=not PM2, deactivate_PP5_BP6=not PP5_BP6, predictor=args.predictor)
        run.vcf_run(args.input, svcf_mode=True)

    else:
        if ".txt" in args.input:
            variant_list = pd.read_csv(args.input, sep="\t", header=None)[0].tolist()
        else:
            variant_list = [args.input]

        run = AAVC(variant_list, deactivate_PM2=not PM2, deactivate_PP5_BP6=not PP5_BP6, clock=True)
        run.run()

    # close DB connection
    connection.close()


