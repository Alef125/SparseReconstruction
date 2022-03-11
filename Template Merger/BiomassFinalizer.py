"""
This code, finalizes the stoichiomety matrix and bounds based on the biomass reaction necessary for growth constraints.
"""

import pandas as pd
import json
import os


class BiomassFinalizer:
    def __init__(self,
                 template_lower_bounds_filepath: str,
                 template_upper_bounds_filepath: str,
                 stoichiometric_data_filepath: str,
                 template_metabolites_filepath: str,
                 existing_reactions_filepath: str,
                 biomass_template_id: str = None,
                 biomass_composition_filepath: str = None,
                 biomass_growth_threshold: float = 1e-6):
        """
        :param template_lower_bounds_filepath: The filepath for all_lower_bounds.csv
        :param template_upper_bounds_filepath: The filepath for all_upper_bounds.csv
        :param stoichiometric_data_filepath: The filepath for stoichiometric_data.json
        :param template_metabolites_filepath: The filepath for template_metabolites.json
        :param existing_reactions_filepath: A string denoting the filepath for existing_rxns.json file.
        :param biomass_template_id: The ID for the biomass reaction, if present in the template
        :param biomass_composition_filepath: The filepath for biomass_composition.csv
        :param biomass_growth_threshold: The minimum biomass production rate for organism's growth
        """
        # #################################################################
        self.all_template_reactions = None
        self.template_lower_bounds_filepath = template_lower_bounds_filepath
        self.template_placed_lower_bounds = None
        self.template_upper_bounds_filepath = template_upper_bounds_filepath
        self.template_placed_upper_bounds = None
        self.load_template_bounds()
        # ##############################################################
        self.stoichiometric_data_filepath = stoichiometric_data_filepath
        self.stoichiometric_data = None
        self.load_stoichiometric_data()
        self.sparse_stoichiometry_matrix = None
        # ##############################################################
        self.template_metabolites_filepath = template_metabolites_filepath
        self.all_template_metabolites = None
        self.load_template_metabolites()
        # ###################################################
        self.existing_reactions_filepath = existing_reactions_filepath
        self.existing_rxns_ids = None
        self.load_existing_reactions()
        # #############################################
        self.biomass_template_id = biomass_template_id
        self.biomass_composition_filepath = biomass_composition_filepath
        self.biomass_growth_threshold = biomass_growth_threshold

    def load_template_bounds(self):
        """
        This method, loads the template bounds .csv file from self.template_lower/upper_bounds_filepath
        :return: -
        """
        self.template_placed_lower_bounds = pd.read_csv(self.template_lower_bounds_filepath)
        self.template_placed_upper_bounds = pd.read_csv(self.template_upper_bounds_filepath)
        self.all_template_reactions = self.template_placed_lower_bounds['ID'].tolist()
        if self.all_template_reactions != self.template_placed_upper_bounds['ID'].tolist():
            print("Your lower and upper bounds are not compatible")
            raise Exception

    def load_stoichiometric_data(self):
        """
        This method, loads the stoichiometric data .json file from self.stoichiometric_data_filepath
        :return: -
        """
        with open(self.stoichiometric_data_filepath, 'r') as json_file:
            self.stoichiometric_data = json.load(json_file)

    def load_template_metabolites(self):
        """
        This method, loads the template_metabolites .json file from self.template_metabolites_filepath
        :return: -
        """
        with open(self.template_metabolites_filepath, 'r') as json_file:
            self.all_template_metabolites = json.load(json_file)

    def load_existing_reactions(self):
        """
        This method, loads the existing reactions .json file from self.existing_reactions_filepath
        :return: -
        """
        with open(self.existing_reactions_filepath, 'r') as json_file:
            self.existing_rxns_ids = json.load(json_file)

    def finalize_biomass_by_composition(self):
        # ToDo (Important one!)
        pass

    def finalize_biomass_by_id(self):
        """
        This method, modified the biomass bounds to be (self.biomass_growth_threshold, 1e6)
        in our growth data of self.template_placed_lower_bounds and self.template_placed_upper_bounds.
        :return: -
        """
        lb_columns = list(self.template_placed_lower_bounds.columns)
        ub_columns = list(self.template_placed_upper_bounds.columns)
        lb_columns.remove('ID')
        ub_columns.remove('ID')
        self.template_placed_lower_bounds.loc[
            self.template_placed_lower_bounds['ID'] == self.biomass_template_id,
            lb_columns
        ] = self.biomass_growth_threshold
        self.template_placed_upper_bounds.loc[
            self.template_placed_upper_bounds['ID'] == self.biomass_template_id,
            ub_columns
        ] = 1e6

    def make_sparse_stoichiometry_matix(self, reactions_index_map: dict, metabolites_index_map: dict):
        """
        This method, parses the self.stoichiometric_data to make the stoichiometry matrix.
        :param reactions_index_map: indexes assigned to the reactions
        :param metabolites_index_map: indexes assigned to the metabolites
        :return: Filling self.sparse_stoichiometry_matrix with a pandas' DataFrame including three columns of
                 met_id, rxn_id, and coeff.
        """
        rows_indexes = []
        cols_indexes = []
        coefficients = []
        for rxn_id, metabolites in self.stoichiometric_data.items():
            for metabolite, stoich_coeff in metabolites.items():
                rows_indexes.append(metabolites_index_map[metabolite])
                cols_indexes.append(reactions_index_map[rxn_id])
                coefficients.append(stoich_coeff)
        self.sparse_stoichiometry_matrix = pd.DataFrame({'met_id': rows_indexes,
                                                         'rxn_id': cols_indexes,
                                                         'coeff': coefficients})

    def save_final_data(self, folder_to_save: str):
        """
        This method, saves all 5 final files in the folder_to_save:
            1. "L.csv": finalized self.template_placed_lower_bounds
            2. "U.csv": finalized self.template_placed_lower_bounds
            3. "reactions_index_map.json": indexes assigned to the reactions
            4. "metabolites_index_map.json": indexes assigned to the metabolites
            5. "S.csv": finalized self.sparse_stoichiometry_matrix
        :param folder_to_save: The folder to save final files.
        :return: -
        """
        if not os.path.exists(folder_to_save):
            os.makedirs(folder_to_save)
        self.template_placed_lower_bounds = self.template_placed_lower_bounds.drop(columns=['ID'])
        self.template_placed_upper_bounds = self.template_placed_upper_bounds.drop(columns=['ID'])
        reactions_index_map = {k: v for v, k in enumerate(self.all_template_reactions)}
        metabolites_index_map = {k: v for v, k in enumerate(self.all_template_metabolites)}
        internal_rxns_indexes = [reactions_index_map[rxn_id] for rxn_id in self.existing_rxns_ids]
        self.make_sparse_stoichiometry_matix(reactions_index_map=reactions_index_map,
                                             metabolites_index_map=metabolites_index_map)
        # ##################################  Saving ####################################
        self.template_placed_lower_bounds.to_csv(folder_to_save + 'L.csv', index=False)
        self.template_placed_upper_bounds.to_csv(folder_to_save + 'U.csv', index=False)
        with open(folder_to_save + 'existing_reactions.json', 'w') as file:
            json.dump(internal_rxns_indexes, file)
        with open(folder_to_save + 'reactions_index_map.json', 'w') as file:
            json.dump(reactions_index_map, file)
        with open(folder_to_save + 'metabolites_index_map.json', 'w') as file:
            json.dump(metabolites_index_map, file)
        self.sparse_stoichiometry_matrix.to_csv(folder_to_save + 'S.csv', index=False)

    def finalize_and_save_data(self, folder_to_save: str):
        """
        This method, chooses the right function to finalize date based on the biomass information.
        Note that composition information has a priority over the biomass id.
        Also, if neither of those parameters are defined, this method raises an Exception.
        :param folder_to_save: The folder to save L.csv, U.csv, and S.csv
        :return: -
        """
        if self.biomass_composition_filepath:
            self.finalize_biomass_by_composition()
        else:
            if self.biomass_template_id:
                self.finalize_biomass_by_id()
            else:
                print("Neither Biomass_composition nor Biomass_id are defined")
                raise Exception
        self.save_final_data(folder_to_save=folder_to_save)


lbs_filepath = "../../Data/Palsson B.Subtilis Reconstruction/Micro-Template Placed Bounds/lower_bounds.csv"
ubs_filepath = "../../Data/Palsson B.Subtilis Reconstruction/Micro-Template Placed Bounds/upper_bounds.csv"
stoich_filepath = "../../Data/Palsson B.Subtilis Reconstruction/Microbial Template/Microbial Stoichiometric Data.json"
mets_filepath = "../../Data/Palsson B.Subtilis Reconstruction/Microbial Template/Microbial Template Metabolites.json"
existing_rxns_filepath = "../Data/Palsson B.Subtilis Reconstruction/existing_rxns.json"
obj = BiomassFinalizer(
    template_lower_bounds_filepath=lbs_filepath,
    template_upper_bounds_filepath=ubs_filepath,
    stoichiometric_data_filepath=stoich_filepath,
    template_metabolites_filepath=mets_filepath,
    existing_reactions_filepath=existing_rxns_filepath,
    biomass_template_id="Growth",  # "BIOMASS_Ec_iAF1260_core_59p81M",
    biomass_growth_threshold=1e-2)

obj.finalize_and_save_data(folder_to_save="../Data/Palsson B.Subtilis Reconstruction/Microbial Final Data/")
