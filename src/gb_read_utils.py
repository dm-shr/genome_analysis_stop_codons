from typing import List, Optional, Dict
from dataclasses import dataclass

import yaml

from Bio import SeqIO
from Bio.Seq import Seq


@dataclass
class Record:
    """
    A class representing a GenBank (GB) record with extracted genome features.

    Attributes:
        tax_id (str): The taxonomic identifier of the organism.
        taxonomy (str): The full taxonomic hierarchy of the organism.
        genes (list, optional): A list containing the names or identifiers of the genes. Default is None.
        sequences (list, optional): A list containing the DNA sequences corresponding to the genes.
            Each entry should correspond to a gene in the 'genes' list. Default is None.
        ends (list, optional): A list containing information about the ends of the sequences or genes.
            Each entry should correspond to a gene in the 'genes' list. Default is None.
        frameshifts (list, optional): A list containing information about about the ends with -1 frameshift
          within the genes.
            Each entry should correspond to a gene in the 'genes' list. Default is None.
        is_outlier (bool, optional): A boolean flag indicating whether this record is considered an outlier in
            Default is False.
    """
    tax_id: str
    entry: str
    taxonomy: str
    is_odd: Optional[bool] = False
    _genes: Optional[List[str]] = None
    _sequences: Optional[List[str]] = None
    _locations: Optional[List[str]] = None
    _ends: Optional[List[str]] = None
    _frameshifts: Optional[List[str]] = None
    _is_outlier: Optional[bool] = False
    

    @staticmethod
    def trim(seq: Seq, frameshift: bool = False) -> str:
        """
        Given an ORF sequence as a Seq object,
        iterate over it and find its last nucleotides:
        could be 3 letters, if ORF length can be divided by 3,
        or 1 or 2 letters.
        """
        length = len(seq)

        if length % 3 == 0:
            index = -3
        elif length % 3 == 1:
            index = -1
        else:
            index = -2

        if frameshift:
            index += -1
        
        return str(seq[index:])

    @staticmethod
    def frameshift(seq: Seq) -> str:
        """
        If last in-frame nucleotides are within a specified set,
        return them.
        """
        if Record.trim(seq) in {'AGA', 'AGG', 'AG', 'A'}:
            return Record.trim(seq, frameshift=True)

        return 'NaN'
        
    def check_outlier(self, canonical_names: Dict[str, str], genes_count: int = 13) -> None:
        """
        Set is_outlier to True if conditions are not met (gene names and genes count).
        Rename gene names to the canonical ones if possible.
        """
        new_names = []
        for gene_name in self._genes:
            try:
                new_names.append(canonical_names[gene_name])
            except KeyError:
                new_names.append(gene_name)
                self._is_outlier = True
        
        self._genes = new_names

        if len(self._genes) != genes_count:
            self._is_outlier = True
    
    def extend(
            self,
            gene_name: str,
            seq: Seq,
            location: str   
    ) -> None:
        "Once a new feature (CDS) is extracted, add it to a class"
        # control for the sizes of the arrays
        check_tuple = (
            self._genes, 
            self._sequences, 
            self._locations,
            self._ends, 
            self._frameshifts
            )
        
        if not all(check_tuple):
            if any(check_tuple):
                raise ValueError(
                    f'tax_id = {self.tax_id}, feature arrays of different length.'\
                    'check size of genes, sequnces, ends, frameshifts'
                )
            
            self._genes = []
            self._sequences = []
            self._locations = []
            self._ends = []
            self._frameshifts = []
            
        # if sequence is not defined (not provided in the genome),
        # let it be empty
        if not seq.defined:
                self.is_odd = True
                seq = ''
        
        # add the elements to the arrays
        # gene name
        self._genes.append(gene_name)

        # full sequence
        self._sequences.append(str(seq))

        # its genome coordinates
        self._locations.append(location)

        # last 1-3 in-frame nucleotides
        self._ends.append(Record.trim(seq))

        # if ORF ends in 'AGA', 'AGG', 'AG', 'A',
        # append -1 frameshift result,
        # else NaN
        self._frameshifts.append(Record.frameshift(seq))
    
    def sort(self):
        "Sort the arrays so that gene names are in aplhabetical order"
        # Zip lists together
        zipped_lists = zip(
            self._genes, 
            self._sequences, 
            self._locations,
            self._ends, 
            self._frameshifts
            )

        # Sort by genes (as it is the first element)
        sorted_lists = sorted(zipped_lists)

        # Unzip lists
        self._genes,\
        self._sequences,\
        self._locations,\
        self._ends,\
        self._frameshifts = zip(*sorted_lists)
    
    def write_tsv(
            self, 
            filepath: str, 
            include_taxonomy: bool = True,
            short_output: bool = True
            ) -> None:
        "writes the object into tsv file"

        if short_output:
            sequences = self._locations
        else:
            sequences = self._sequences
        
        taxonomy = 'NaN'
        if include_taxonomy:
            taxonomy = self.taxonomy

        with open(filepath, 'a') as file:
            iter_object = zip(
                self._genes, 
                sequences, 
                self._ends, 
                self._frameshifts
                )
            for gene_name, sequence, end, frameshift in iter_object:
                string_to_write = "\t".join(
                    (
                        self.tax_id,
                        self.entry,
                        taxonomy,
                        str(self._is_outlier),
                        str(self.is_odd),
                        gene_name,
                        sequence,
                        end,
                        frameshift
                    ))
                file.write(string_to_write + '\n')



def write_genome_summary(config_path: str = './config.yml'):
    """given gb files, write a tsv table with the columns for each:
    tax_id,
    entry,
    taxonomy,
    is_outlier,
    is_odd,
    gene_name,
    sequence,
    end,
    frameshift.
    """
    # Load the config and initialize the params
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)

    # for a final tsv table
    output_path = config['output_path']

    # input gb file path
    gb_file_path = config['gb_file_path']

    # allowed taxonomy path to filter out target species
    taxonomy_path = config['taxonomy_path']
    # allowed names on mitochondrial genes
    genes_names_dict = {
        alias.upper(): name
        for alias, name in config['genes_names_dict'].items()
    }

    # initiliaze genbank object
    gb_dict = SeqIO.index(gb_file_path, 'genbank')

    # gb_dicts = [SeqIO.index(gb_file, 'genbank') for gb_file in gb_files]

    # # iterate over GB files:
    # for gb_dict in gb_dicts:
    #     # break

    # iterate over GB records
    for entry, data in gb_dict.items():

        # odd indicates non-provided genome sequences,
        # or possibly wrong taxonomy annotation
        is_odd = False

        # filter out non-relevant species
        if data.annotations['taxonomy'][:len(taxonomy_path)] != taxonomy_path:
            # if all(
            #     [
            #         taxon.lower() != 'vertebrata' for taxon in data.annotations['taxonomy']
            #         ]
            #         ):
            #     continue

            # is_odd = True
            continue
            
        tax_id = 'NaN'
        # try to infer taxonomy id from the record data
        if 'db_xref' in data.features[0].qualifiers:
            tax_id = data.features[0].qualifiers['db_xref'][0].split(':')[1]

        # create Record object
        taxonomy = ':'.join(data.annotations['taxonomy'])
        record = Record(tax_id, entry, taxonomy, is_odd=is_odd)

        # iterate over features
        for feature in data.features[1:]:

            # filter out non-CDS features
            if feature.type != 'CDS':
                continue
            
            # get gene sequence
            seq = feature.extract(data.seq)

            # get sequence location on the genome coordinates
            location = str(feature.location)

            # get gene name
            if 'gene' in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0].upper()
            else:
                gene_name = feature.qualifiers['product'][0].upper()

            record.extend(gene_name=gene_name, seq=seq, location=location)

        # sort the arrays by gene_name
        record.sort()
        # check if the record is an outlier 
        # (not equal to 13 CDS, non-canonical names)
        # and renames the genes to the canonical names
        record.check_outlier(canonical_names=genes_names_dict)
        # write the record to a file
        record.write_tsv(output_path, include_taxonomy=False)


if __name__ == '__main__':
    write_genome_summary()
