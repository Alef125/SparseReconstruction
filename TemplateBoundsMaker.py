"""
TemplateBoundsMaker
This code, provides a class to place organism's bounds into the template.
"""

import pandas as pd
import warnings
import os
import json


def make_cell_to_cell_translation_dict(key_cell: str, value_cell: str) -> dict:
    """
    :param key_cell: The string content of a cell in the template_bounds[input_reactions_nomenclature],
                     like "[rxn1; rxn2; rxn3]"
    :param value_cell: The string content of a cell in the template_bounds[template_reactions_nomenclature],
                       like "[rxn4; rxn5; rxn6]"
    :return: translation_dict: Dict of all possible translations, like:
                        {rxn1: [rxn4, rxn5, rxn6], rxn2: [rxn4, rxn5, rxn6], rxn3: [rxn4, rxn5, rxn6]}
    """
    translation_dict = {}
    if pd.isna(key_cell):
        return {}
    key_reactions = key_cell.split("; ")
    if pd.isna(value_cell):
        value_reactions = []
    else:
        value_reactions = value_cell.split("; ")
    for key_reaction in key_reactions:
        translation_dict[key_reaction] = value_reactions
    return translation_dict


class TemplateBoundsMaker:
    def __init__(self,  # ToDo: Biomass!!!
                 lower_bounds_filepaths: list,
                 upper_bounds_filepaths: list,
                 internal_rxns_filepath: str,
                 template_bounds_filepath: str,
                 input_reactions_nomenclature: str = None,
                 template_reactions_nomenclature: str = None,
                 reactions_translation_filepath: str = None):
        """
        :param lower_bounds_filepaths: List of paths for all .csv lower_bounds to be merged and placed on the template.
        :param upper_bounds_filepaths: List of paths for all .csv upper_bounds to be merged and placed on the template.
        :param internal_rxns_filepath: A string denoting the filepath for internal_rxns_bounds.csv file.
        :param template_bounds_filepath: The path for the template_bounds.csv file
        :param reactions_translation_filepath: The path for the reactions_translation.csv file
        :param input_reactions_nomenclature: The column name in the translation_file
                                             corresponding to the reaction names in the lower/upper_bounds files
        :param template_reactions_nomenclature: The column name in the translation_file
                                                corresponding to the template_bounds file
        """
        self.total_reactions_list = None
        self.lower_bounds_filepaths = lower_bounds_filepaths
        self.lower_bounds_df = None
        self.read_and_merge_lower_bounds()
        self.upper_bounds_filepaths = upper_bounds_filepaths
        self.upper_bounds_df = None
        self.read_and_merge_upper_bounds()
        # ###################################################
        self.internal_rxns_filepath = internal_rxns_filepath
        self.internal_rxns_ids = None
        self.load_internal_reactions()
        # ######################################################
        self.template_bounds_filepath = template_bounds_filepath
        self.template_bounds = None
        self.all_template_reactions = None
        self.load_template_bounds()
        # #############################
        self.input_reactions_nomenclature = input_reactions_nomenclature
        self.template_reactions_nomenclature = template_reactions_nomenclature
        # #############################
        self.reactions_translation_filepath = reactions_translation_filepath
        self.translation_dict = None
        self.make_translation_dict()
        # ###############################
        self.internal_rxns_temp = []
        self.template_placed_lower_bounds = None
        self.template_placed_upper_bounds = None

    def read_and_merge_lower_bounds(self):  # ToDo: columns names (and confidence)
        """
        This method, reads all lower_bound files of self.lower_bounds_filepaths and merges them together
        :return: filling the self.lower_bounds_df
        """
        do_initiate = True
        for lower_bounds_filepath in self.lower_bounds_filepaths:
            lower_bounds_file = pd.read_csv(lower_bounds_filepath)
            if do_initiate:
                self.total_reactions_list = lower_bounds_file['ID'].tolist()
                self.lower_bounds_df = lower_bounds_file
                do_initiate = False
            else:
                lower_bounds_file = lower_bounds_file.drop('ID', 1)  # 1: column, 0:row
                # Default inner join:
                self.lower_bounds_df = pd.merge(self.lower_bounds_df, lower_bounds_file,
                                                left_index=True, right_index=True)

    def read_and_merge_upper_bounds(self):
        """
        This method, reads all upper_bound files of self.upper_bounds_filepaths and merges them together
        :return: filling the self.upper_bounds_df
        """
        do_initiate = True
        for upper_bounds_filepath in self.upper_bounds_filepaths:
            upper_bounds_file = pd.read_csv(upper_bounds_filepath)
            if do_initiate:
                self.total_reactions_list = upper_bounds_file['ID'].tolist()
                self.upper_bounds_df = upper_bounds_file
                do_initiate = False
            else:
                upper_bounds_file = upper_bounds_file.drop('ID', 1)  # 1: column, 0:row
                # Default inner join:
                self.upper_bounds_df = pd.merge(self.upper_bounds_df, upper_bounds_file,
                                                left_index=True, right_index=True)

    def load_internal_reactions(self):
        """
        This method, loads the internal reactions bounds .csv file from self.internal_rxns_filepath
        :return: -
        """
        internal_rxns_df = pd.read_csv(self.internal_rxns_filepath)
        self.internal_rxns_ids = internal_rxns_df['ID'].tolist()

    def load_template_bounds(self):
        """
        This method, loads the template bounds .csv file from self.template_bounds_filepath
        :return: -
        """
        self.template_bounds = pd.read_csv(self.template_bounds_filepath)
        self.all_template_reactions = self.template_bounds['ID'].tolist()

    def make_translation_dict(self):
        """
        This method, reads reactions_translation_file.csv from self.reactions_translation_filepath,
            then makes the self.translation_dict dictionary with keys as reactions names in
            self.input_reactions_nomenclature and values as reactions names in self.template_reactions_nomenclature
        """
        if self.reactions_translation_filepath:
            self.translation_dict = {}
            translation_file = pd.read_csv(self.reactions_translation_filepath)
            if self.input_reactions_nomenclature:
                if self.input_reactions_nomenclature not in translation_file.columns:
                    print("input_reactions_nomenclature does not exist in translation_file")  # ToDo: to Exception
                    raise Exception
            if self.template_reactions_nomenclature:
                if self.template_reactions_nomenclature not in translation_file.columns:
                    print("template_reactions_nomenclature does not exist in translation_file")
                    raise Exception
            for index, row in translation_file.iterrows():
                key_cell = row[self.input_reactions_nomenclature]
                value_cell = row[self.template_reactions_nomenclature]
                translation_dict = make_cell_to_cell_translation_dict(key_cell=key_cell, value_cell=value_cell)
                self.translation_dict.update(translation_dict)
                translation_dict2 = make_cell_to_cell_translation_dict(key_cell=value_cell, value_cell=value_cell)
                self.translation_dict.update(translation_dict2)  # ToDo (important): Dirty addition

    def translate_internal_reactions(self):
        """
        This method, finds corresponding ids for self.internal_rxns_ids in the template.
        :return: Filling self.internal_rxns_temp.
        """
        for rxn_base_id in self.internal_rxns_ids:
            if rxn_base_id in self.translation_dict.keys():
                rxn_temp_ids = self.translation_dict[rxn_base_id]
                if not rxn_temp_ids:
                    warn_text = "Bad translation for " + rxn_base_id  # 161
                    warnings.warn(warn_text)
                else:
                    if rxn_base_id in rxn_temp_ids:  # 745
                        self.internal_rxns_temp.append(rxn_base_id)
                    else:  # Happens for 8 reactions, ToDo: Finding Best id is so rudimentary
                        self.internal_rxns_temp.append(rxn_temp_ids[0])
            else:
                warn_text = "The internal reaction " + rxn_base_id + " is not included in your translation file."  # 26
                warnings.warn(warn_text)

    def make_template_lower_bounds(self):
        """
        This method, overrides the template lower bounds by existing organism's bounds based on self.translation_dict.
        :return: Filling self.template_placed_lower_bounds.
        """
        # ########## Initiate self.template_placed_lower_bounds with default Bounds #############
        all_data_columns = list(self.lower_bounds_df.columns)
        all_data_columns.remove('ID')
        self.template_placed_lower_bounds = pd.DataFrame({'ID': self.all_template_reactions})
        for data_column in all_data_columns:
            self.template_placed_lower_bounds[data_column] = self.template_bounds['Lower Bound']
        # ########### Finding corresponding rows and modify them in template rows ##############
        for index, row in self.lower_bounds_df.iterrows():
            if row['ID'] in self.translation_dict.keys():
                template_ids = self.translation_dict[row['ID']]
                in_template = False
                for template_id in template_ids:
                    if template_id in self.all_template_reactions:
                        self.template_placed_lower_bounds.loc[
                            self.template_placed_lower_bounds['ID'] == template_id,
                            all_data_columns
                        ] = row[all_data_columns].tolist()
                        in_template = True
                if not in_template:
                    warning_text = "The reaction with ID " + row['ID'] + \
                                   " does not have any corresponding reaction in the template"
                    warnings.warn(warning_text)  # 16 warnings (148)
            else:
                warn_text = "The reaction " + row['ID'] + " is not included in your translation file."
                warnings.warn(warn_text)  # 451 warnings (902)(251)

    def make_template_upper_bounds(self):
        """
        This method, overrides the template upper bounds by existing organism's bounds based on self.translation_dict.
        :return: Filling self.template_placed_upper_bounds.
        """
        # ########## Initiate self.template_placed_upper_bounds with default Bounds #############
        all_data_columns = list(self.upper_bounds_df.columns)
        all_data_columns.remove('ID')
        self.template_placed_upper_bounds = pd.DataFrame({'ID': self.all_template_reactions})
        for data_column in all_data_columns:
            self.template_placed_upper_bounds[data_column] = self.template_bounds['Upper Bound']
        # ########### Finding corresponding rows and modify them in template rows ##############
        for index, row in self.upper_bounds_df.iterrows():
            if row['ID'] in self.translation_dict.keys():
                template_ids = self.translation_dict[row['ID']]
                in_template = False
                for template_id in template_ids:
                    if template_id in self.all_template_reactions:
                        self.template_placed_upper_bounds.loc[
                            self.template_placed_upper_bounds['ID'] == template_id,
                            all_data_columns
                        ] = row[all_data_columns].tolist()
                        in_template = True
                if not in_template:
                    warning_text = "The reaction with ID " + row['ID'] + \
                                   " does not have any corresponding reaction in the template"
                    warnings.warn(warning_text)  # 16 warnings
            else:
                warn_text = "The reaction " + row['ID'] + " is not included in your translation file."
                warnings.warn(warn_text)  # 451 warnings (251)

    def save_existing_reaction(self, path_to_save: str):
        """
        :param path_to_save: File path to save self.internal_rxns_temp
        :return: -
        """
        with open(path_to_save, 'w') as file:
            json.dump(self.internal_rxns_temp, file)

    def save_final_bounds(self, folder_to_save: str):
        """
        :param folder_to_save: Folder path to save final lower and upper bound
        :return: -
        """
        if not os.path.exists(folder_to_save):
            os.makedirs(folder_to_save)
        self.template_placed_lower_bounds.to_csv(folder_to_save + 'lower_bounds.csv', index=False)
        self.template_placed_upper_bounds.to_csv(folder_to_save + 'upper_bounds.csv', index=False)


g_lb_filepaths = ["../Data/Palsson B.Subtilis Reconstruction/Util Bounds/g_lower_bounds.csv",
                  "../Data/Palsson B.Subtilis Reconstruction/KO Bounds/g_lower_bounds.csv"]
g_ub_filepaths = ["../Data/Palsson B.Subtilis Reconstruction/Util Bounds/g_upper_bounds.csv",
                  "../Data/Palsson B.Subtilis Reconstruction/KO Bounds/g_upper_bounds.csv"]

ng_lb_filepaths = ["../Data/Palsson B.Subtilis Reconstruction/Util Bounds/ng_lower_bounds.csv",
                   "../Data/Palsson B.Subtilis Reconstruction/KO Bounds/ng_lower_bounds.csv"]
ng_ub_filepaths = ["../Data/Palsson B.Subtilis Reconstruction/Util Bounds/ng_upper_bounds.csv",
                   "../Data/Palsson B.Subtilis Reconstruction/KO Bounds/ng_upper_bounds.csv"]

micro_template_filepath = "../Data/Palsson B.Subtilis Reconstruction/Microbial Template/Microbial Universal Bounds.csv"
g_obj = TemplateBoundsMaker(
    lower_bounds_filepaths=g_lb_filepaths,
    upper_bounds_filepaths=g_ub_filepaths,
    internal_rxns_filepath="../Data/Palsson B.Subtilis Reconstruction/Internal_Rxns_Bounds.csv",
    template_bounds_filepath=micro_template_filepath,
    input_reactions_nomenclature="KBase Abbr",
    template_reactions_nomenclature="BiGG id",
    reactions_translation_filepath="../Data/Palsson B.Subtilis Reconstruction/Reactions Translation.csv")
g_obj.translate_internal_reactions()
g_obj.make_template_lower_bounds()
g_obj.make_template_upper_bounds()
g_obj.save_existing_reaction(path_to_save="../Data/Palsson B.Subtilis Reconstruction/existing_rxns.json")
g_obj.save_final_bounds(folder_to_save="../Data/Palsson B.Subtilis Reconstruction/Micro-Template Placed Bounds/")

# ng_obj = TemplateBoundsMaker(
#     lower_bounds_filepaths=ng_lb_filepaths,
#     upper_bounds_filepaths=ng_ub_filepaths,
#     template_bounds_filepath="../Data/Palsson B.Subtilis Reconstruction/Universal Bounds.csv",
#     input_reactions_nomenclature="BiGG id",
#     template_reactions_nomenclature="BiGG id",
#     reactions_translation_filepath="../Data/Palsson B.Subtilis Reconstruction/Reactions Translation.csv")
# ng_obj.make_template_lower_bounds()
# ng_obj.make_template_upper_bounds()
# ng_obj.save_final_bounds(folder_to_save="../Data/Palsson B.Subtilis Reconstruction/Template Placed Bounds/")
