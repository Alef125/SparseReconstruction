"""
ReactionsTranslation
The main class in this script, Translator, builds an object to translate ids of reactions from an input_nomenclature
into a template_nomenclature.
"""

import warnings
import pandas as pd


def get_cell_str_to_list(cell_str: str):
    """
    This function, parses the string content of a template_bounds cell and into a list of ids.
    The overall format of the string upon brackets ([]), quotations('', ""), or delimiter (, or ;)
    is recognized Automatically.
    Acceptable formats: ['id1', 'id2', 'id3'] - ['id1'; 'id2'; 'id3']
                        ["id1", "id2", "id3"] - ["id1"; "id2"; "id3"]
                        'id1', 'id2', 'id3' - 'id1'; 'id2'; 'id3'
                        "id1", "id2", "id3" - "id1"; "id2"; "id3"
    :param cell_str: The string content of a template_bounds cell
    :return: ids_list
    """
    if cell_str[0] == '[' and cell_str[-1] == ']':
        cell_str = cell_str[1:-1]
    if ',' in cell_str:
        delimiter = "; "
    else:
        delimiter = "; "
    ids_raw_list = cell_str.split(delimiter)
    ids_list = []
    for raw_id in ids_raw_list:
        if raw_id[0] == '\'' and raw_id[-1] == '\'':
            ids_list.append(raw_id[1:-1])
        elif raw_id[0] == '\"' and raw_id[-1] == '\"':
            ids_list.append(raw_id[1:-1])
        else:
            ids_list.append(raw_id)
    return ids_list


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
    key_reactions = get_cell_str_to_list(cell_str=key_cell)
    if pd.isna(value_cell):
        value_reactions = []
    else:
        value_reactions = get_cell_str_to_list(cell_str=value_cell)
    for key_reaction in key_reactions:
        translation_dict[key_reaction] = value_reactions
    return translation_dict


def make_general_translation_dict(reactions_translation_filepath: str,
                                  input_reactions_nomenclature: str,
                                  template_reactions_nomenclature: str) -> dict:
    """
    This method, reads reactions_translation_file.csv from reactions_translation_filepath,
        then makes the general_translation_dict with keys as reactions names in input_reactions_nomenclature
        and values as reactions names in template_reactions_nomenclature
    :param reactions_translation_filepath: Filepath for the reactions_translation_file.csv
    :param input_reactions_nomenclature: Source IDs for translation
    :param template_reactions_nomenclature: Destination IDs for translation
    :return: general_translation_dict
    """
    translation_file = pd.read_csv(reactions_translation_filepath)
    if input_reactions_nomenclature:
        if input_reactions_nomenclature not in translation_file.columns:
            raise KeyError("input_reactions_nomenclature does not exist in translation_file")
    if template_reactions_nomenclature:
        if template_reactions_nomenclature not in translation_file.columns:
            raise KeyError("template_reactions_nomenclature does not exist in translation_file")
    general_translation_dict = {}
    for index, row in translation_file.iterrows():
        key_cell = row[input_reactions_nomenclature]
        value_cell = row[template_reactions_nomenclature]
        translation_dict = make_cell_to_cell_translation_dict(key_cell=key_cell, value_cell=value_cell)
        general_translation_dict.update(translation_dict)
        translation_dict2 = make_cell_to_cell_translation_dict(key_cell=value_cell, value_cell=value_cell)
        general_translation_dict.update(translation_dict2)  # ToDo (important): Dirty addition
    return general_translation_dict


class Translator:
    def __init__(self,
                 reactions_translation_filepath: str,
                 input_reactions_nomenclature: str,
                 template_reactions_nomenclature: str,
                 list_of_input_reactions: list,
                 list_of_template_reactions: list):
        """
        :param reactions_translation_filepath: The path for the reactions_translation.csv file
        :param input_reactions_nomenclature: The column name in the translation_file
                                             corresponding to the reaction names in the lower/upper_bounds files
        :param template_reactions_nomenclature: The column name in the translation_file
                                                corresponding to the template_bounds file
        :param list_of_input_reactions: List of desired reactions to be included in the translation (as keys)
        :param list_of_template_reactions: List of valid ids as translation outputs
        """
        # ##############################################################
        self.input_reactions_nomenclature = input_reactions_nomenclature
        self.template_reactions_nomenclature = template_reactions_nomenclature
        # #######################################################
        self.list_of_input_reactions = list_of_input_reactions
        self.list_of_template_reactions = list_of_template_reactions
        # ################################################################
        self.reactions_translation_filepath = reactions_translation_filepath
        self.translation_dict = {}
        self.make_translation_dict()

    def make_translation_dict(self):
        """
        This method, finds corresponding reactions of source_reactions in the destination_reactions.
        :return: Filling self.translation_dict
        """
        # ############ Building general_translation_dict ############
        general_translation_dict = {}
        if self.reactions_translation_filepath:
            general_translation_dict = make_general_translation_dict(
                reactions_translation_filepath=self.reactions_translation_filepath,
                input_reactions_nomenclature=self.input_reactions_nomenclature,
                template_reactions_nomenclature=self.template_reactions_nomenclature
            )
        # ######## Filling self.translation_dict by self.total_reactions_list one by one #########
        for base_reaction_id in self.list_of_input_reactions:
            if base_reaction_id in self.list_of_template_reactions:
                # No need for translation
                self.translation_dict[base_reaction_id] = [base_reaction_id]
            else:
                # Translation needed
                if base_reaction_id in general_translation_dict.keys():
                    # Reaction exists in the general_translation_dict
                    translated_ids = general_translation_dict[base_reaction_id]
                    template_ids = [_rxn_id for _rxn_id in translated_ids
                                    if _rxn_id in self.list_of_template_reactions]
                    if template_ids:
                        # Intersection of translated_ids and self.all_template_reactions is not empty
                        self.translation_dict[base_reaction_id] = template_ids
                    else:
                        # The intersection of translated_ids and self.all_template_reactions is empty (# = 172 of 1165)
                        warning_text = "The reaction with ID " + base_reaction_id + \
                                       " does not have any corresponding reaction in the template"
                        warnings.warn(warning_text)
                        self.translation_dict[base_reaction_id] = []
                else:
                    # Reaction doesn't exist in the general_translation_dict (# = 43 of 1165)
                    warn_text = "The reaction " + base_reaction_id + " is not included in your translation file."
                    warnings.warn(warn_text)
                    self.translation_dict[base_reaction_id] = []
                """
                Note: Exchange reactions map have been manually added to the KBase_Translation.csv file 
                      from row 43776 to 43973. 
                      For more accurate changes, you should modify this end rows or modify the 
                      make_general_translation_dict method to accept an additional exchange_mapping file.
                """

    def translate(self, input_id: str) -> list:
        """
        This method, translates desirable reaction ids into the template_reactions_nomenclature
        :param input_id: A reaction id in the input_reactions_nomenclature to be translated
        :return: Translated ids
        """
        return self.translation_dict[input_id]
