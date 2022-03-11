"""
ReactionsKOMaker.py
This code, translates the genes knocking-out into reactions shutting off, using organism_GPR.json,
and integrates it with the growth information in a specific medium saved in the GenesKO.json file
"""

import json


class ReactionsKOMaker:
    def __init__(self, organism_gpr_filepath, genes_ko_growth_filepath):
        self.organism_gpr_filepath = organism_gpr_filepath
        self.organism_gpr = None
        self.read_organism_gpr()
        self.genes_ko_growth_filepath = genes_ko_growth_filepath
        self.genes_ko_growth_list = None
        self.read_genes_ko_growth_file()

    def read_organism_gpr(self):
        """
        This method, reads organism_gpr (dict) from self.organism_gpr_filepath
        """
        with open(self.organism_gpr_filepath, 'r') as json_file:
            self.organism_gpr = json.load(json_file)

    def read_genes_ko_growth_file(self):
        """
        This method, reads genes_ko_growth_list (list of dicts) from self.genes_ko_growth_filepath
        """
        with open(self.genes_ko_growth_filepath, 'r') as json_file:
            self.genes_ko_growth_list = json.load(json_file)

    def make_reactions_ko_growth(self, filepath_to_save):
        """
        :param filepath_to_save: The path to save the .json file
        :return: Nothing, Saves the list of dictionaries as a .json file
        """
        reactions_ko_dicts = []
        for ko_growth in self.genes_ko_growth_list:
            ko_gene = ko_growth['ko_gene_id']
            if ko_gene in self.organism_gpr.keys():
                shut_reactions = self.organism_gpr[ko_gene]
                reactions_ko_dicts.append(
                    {'ko_rxns_ids': shut_reactions,
                     'medium': ko_growth['medium'],
                     'growth': ko_growth['growth']}
                )
            # else, this gene does not shut off any reaction (e.g. in GPAOr) or is not specified in any complex
        with open(filepath_to_save, 'w', encoding='utf-8') as f:
            json.dump(reactions_ko_dicts, f, ensure_ascii=False, indent=4)


obj = ReactionsKOMaker(organism_gpr_filepath="../../Data/Palsson B.Subtilis Reconstruction/Organism GPR.json",
                       genes_ko_growth_filepath="../../Data/Palsson B.Subtilis Reconstruction/Genes KO Growth.json")
obj.make_reactions_ko_growth(filepath_to_save="../../Data/Palsson B.Subtilis Reconstruction/Reactions KO Growth.json")
