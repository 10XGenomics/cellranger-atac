"""
This implements reference "metadata" operations for single species and barnyard references.
The names for the contigs in the reference are formated as either: "chr<XXX>" or "<species>_chr<XXX>".
This operates on a file, "contigs.json" found in the reference fasta directory that describes
the layout.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

import json
import os.path


class ReferenceManager:
    def __init__(self, path):
        self.path = path

        # Load reference metadata
        if self.metadata_file:
            with open(self.metadata_file, 'r') as infile:
                self.metadata = json.load(infile)
        else:
            self.metadata = None

        # Load contigs defs
        # Check contigs is a valid JSON first before loading
        try:
            with open(self.contig_definitions, 'r') as infile:
                self.contigs = json.load(infile)
                if 'non_nuclear_contigs' in self.contigs and self.contigs['non_nuclear_contigs'] is None:
                    self.contigs['non_nuclear_contigs'] = []

        except ValueError, e:
            raise Exception("Malformed json in '%s': %s" % (self.contig_definitions, e))

        self.sex_chromosomes = {}
        # specification of sex_chromosomes is optional
        if "sex_chromosomes" in self.contigs:
            if not isinstance(self.contigs["sex_chromosomes"], dict):
                raise Exception("Malformed entry for sex_chromosomes in {}".format(self.contig_definitions))
            for species in self.contigs["sex_chromosomes"]:
                for contig in self.contigs["sex_chromosomes"][species]:
                    self.sex_chromosomes[contig] = True

        self.all_contigs = []
        self.contig_lengths = {}
        with open(self.fasta_index, 'r') as infile:
            for line in infile:
                line = line.strip().split()
                name, length = line[0], line[1]
                self.all_contigs.append(name)
                self.contig_lengths[name] = int(length)

    def _filepath_or_error(self, subpaths, err_msg):
        filepath = os.path.join(self.path, *subpaths)
        if not os.path.exists(filepath):
            raise IOError(err_msg)
        return filepath

    def _filepath_or_none(self, subpaths):
        filepath = os.path.join(self.path, *subpaths)
        return filepath if os.path.exists(filepath) else None

    @property
    def metadata_file(self):
        return self._filepath_or_none(["metadata.json"])

    @property
    def genome(self):
        filepath = self._filepath_or_error(["genome"],
                                           "Reference is missing genome name file")
        with open(filepath) as infile:
            return infile.readline().strip()

    @property
    def contig_definitions(self):
        return self._filepath_or_error(["fasta", "contig-defs.json"],
                                       "Reference is missing contig definition file")

    @property
    def fasta(self):
        return self._filepath_or_error(["fasta", "genome.fa"],
                                       "Reference is missing genome fasta sequence")

    @property
    def fasta_index(self):
        return self._filepath_or_error(["fasta", "genome.fa.fai"],
                                       "Reference is missing fasta index file")

    @property
    def genes(self):
        return self._filepath_or_error(["genes", "genes.gtf"],
                                       "Reference is missing genes gtf file")

    @property
    def motifs(self):
        return self._filepath_or_none(["regions", "motifs.pfm"])

    @property
    def blacklist_track(self):
        return self._filepath_or_none(["regions", "blacklist.bed"])

    @property
    def ctcf_track(self):
        return self._filepath_or_none(["regions", "ctcf.bed"])

    @property
    def dnase_track(self):
        return self._filepath_or_none(["regions", "dnase.bed"])

    @property
    def enhancer_track(self):
        return self._filepath_or_none(["regions", "enhancer.bed"])

    @property
    def promoter_track(self):
        return self._filepath_or_none(["regions", "promoter.bed"])

    @property
    def tss_track(self):
        return self._filepath_or_none(["regions", "tss.bed"])

    @property
    def transcripts_track(self):
        return self._filepath_or_none(["regions", "transcripts.bed"])

    # TODO: make all these functions properties (and update usage)
    def non_nuclear_contigs(self):
        """ Lists non-nuclear contigs (e.g. mitochondria and chloroplasts) """
        key = "non_nuclear_contigs"
        return self.contigs[key] if key in self.contigs else []

    def get_contig_lengths(self):
        return self.contig_lengths

    def get_vmem_est(self):
        total_contig_len = sum(self.contig_lengths.values())
        return 1.2 * total_contig_len / 1e9

    def primary_contigs(self, species=None, allow_sex_chromosomes=False):
        """List all primary contigs, optionally filter by species."""
        if species is not None:
            if species not in self.contigs["species_prefixes"]:
                raise ValueError("Unknown species")
        return [x for x in self.contigs["primary_contigs"]
                if self.is_primary_contig(x, species=species, allow_sex_chromosomes=allow_sex_chromosomes)]

    def list_species(self):
        return self.contigs["species_prefixes"]

    def is_primary_contig(self, contig_name, species=None, allow_sex_chromosomes=False):
        if contig_name not in self.contigs["primary_contigs"]:
            return False
        if species is not None and self.species_from_contig(contig_name) != species:
            return False
        if allow_sex_chromosomes is False:
            if contig_name in self.sex_chromosomes:
                return False
        return True

    def species_from_contig(self, contig):
        """Given the name of a contig, extract the species."""
        if self.contigs["species_prefixes"] == [""]:
            return ""

        s = contig.split("_")
        if len(s) > 1:
            return s[0]
        else:
            return ""

    def verify_contig_defs(self):
        """Verify that the contig defs are correctly formatted"""
        def check_aos(a):
            if not isinstance(a, list):
                return False

            for k in a:
                if not isinstance(k, unicode):
                    return False
            return True

        # Step 1 does the file exist? -> duplicated functionality in contig_definitions
        if not os.path.exists(self.contig_definitions):
            return "Contig definitions file '%s' is missing" % (self.contig_definitions)

        # Step 2 is it valid JSON? -> implemented before loading during ReferenceManager initiation
        contigs = self.contigs

        # Step 3: species_prefixes must be an array of strings
        species_prefixes = contigs.get("species_prefixes", [u""])
        if not check_aos(species_prefixes):
            return "Species_prefixes must be an array of strings"
        if len(species_prefixes) == 1 and species_prefixes[0] != "":
            return "Cannot specify species prefix for single species references"

        # Step 4: primary contigs must be an array of strings
        if not check_aos(contigs.get("primary_contigs")):
            return "Primary_contigs must be an array of strings"

        # Step 5: prefix contigs can not be prefixes of themselves, and shouldn't contain underscores
        for p1 in species_prefixes:
            for p2 in species_prefixes:
                if p1 != p2 and p1.startswith(p2):
                    return "Species_prefixes are ambiguous. No prefix may be a prefix of another prefix."
        for p in species_prefixes:
            if "_" in p:
                return "species prefix {} contains and underscore. Underscores are not allowed in species prefixes".format(p)

        # step 6: every primary contig and non_nuclear_contigs must be prefixed by a species prefix. Also, every
        # species prefix must be used
        used_prefix = {p: False for p in species_prefixes}
        for c in contigs["primary_contigs"] + contigs["non_nuclear_contigs"]:
            ok = False
            c = str(c)
            for p in species_prefixes:
                if c.startswith(p):
                    used_prefix[p] = True
                    ok = True
                if c == p:
                    return "Species prefix {} matches a contig completely. Species should be a partial prefix to contig".format(p)
            if not ok:
                return "Each primary contig must be prefixed by a species prefix"
        not_used_prefix = [p for p in species_prefixes if not used_prefix[p]]
        if len(not_used_prefix) > 0:
            return "Species prefixes {} are not prefixes for any contigs in the genome fasta".format(",".join(used_prefix))

        # Step 7: if sex_chromosomes are specified
        # they must be a map of maps; Each sub-maps keys must be a primary contig; each sub-maps values must be an integer
        if "sex_chromosomes" in contigs:
            if contigs["sex_chromosomes"] is None or not isinstance(contigs["sex_chromosomes"], dict):
                return "Sex chromosomes must be an object of objects of integers."

            for k in contigs["sex_chromosomes"]:
                v = contigs["sex_chromosomes"][k]
                if not isinstance(v, dict):
                    return "Sex chromosomes must be an object of objects of integers. {} is not an object.".format(k)
                for sx in v:
                    if sx not in contigs["primary_contigs"]:
                        return "Sex chromosomes must be primary contigs. %s is not a primary contig." % sx
                    if not isinstance(v[sx], int):
                        return "Sex chromosome ploidies must be integers."

        # Step 8: every primary contig and non_nuclear_contigs must exist in the reference
        all_fasta_contigs = []

        for line in open(self.fasta_index):
            fields = line.strip().split()
            all_fasta_contigs.append(fields[0])

        for c in contigs["primary_contigs"] + contigs["non_nuclear_contigs"]:
            if c not in (all_fasta_contigs):
                return "Contig %s is in the config definition but not in the genome fasta. Please check if the reference is properly built." % (c)

        # Step 9: there must be a primary contig
        if len(contigs["primary_contigs"]) == 0:
            return "At least one contig must be primary."
        if len(contigs["primary_contigs"]) > 100:
            return "There can be at most 100 primary contigs."

        return None
