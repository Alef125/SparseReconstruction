"""
ActiveNetworkVisualizer

"""

import pandas as pd
import json
from cobra import Model, Metabolite, Reaction
# import networkx as nx
# import matplotlib.pyplot as plt
import cobra
import escher


class ActiveNetworkVisualizer:
    def __init__(self,
                 list_of_active_reactions: list,
                 stoichiometry_dict_filepath: str):
        """
        :param list_of_active_reactions: List of ids for active reactions in a feasible flux
        :param stoichiometry_dict_filepath: The filepath for stoichiometry_dictionary
        """
        # ######################################################
        # ToDo: read list_of_active_reactions from a filepath
        self.list_of_active_reactions = list_of_active_reactions
        # ######################################################
        self.stoichiometry_dict_filepath = stoichiometry_dict_filepath
        self.stoichiometric_dict = None
        self.load_stoichiometric_data()
        # ################################
        self.active_network_stoich = None
        self.active_reactions = []
        self.active_metabolites = []
        self.make_active_network_stoich()
        # ################################
        self.cobra_model = None
        # ######################

    def load_stoichiometric_data(self):
        """
        This method, loads the stoichiometric data .json file from self.stoich_dict_filepath
        :return: -
        """
        with open(self.stoichiometry_dict_filepath, 'r') as json_file:
            self.stoichiometric_dict = json.load(json_file)

    def make_active_network_stoich(self):
        """
        This method, extracts the active portion of the self.stoichiometric_dict
        :return:
        """
        self.stoichiometric_dict["AGPATr_BS"] = {"1ddecg3p_c": -1.0, "ddcaACP_c": -1.0, "ACP_c": 1.0, "pa120_c": 1.0}
        self.stoichiometric_dict["BIOMASS_BS_10"] = self.stoichiometric_dict["Growth"]
        self.stoichiometric_dict["ATPS4r"] = self.stoichiometric_dict["ATPS4rpp"]
        self.stoichiometric_dict["CYOR3m"] = self.stoichiometric_dict["CYO4pp"]
        self.stoichiometric_dict["CYOO3"] = self.stoichiometric_dict["CYO1b"]
        self.stoichiometric_dict["CDGPT_BS"] = self.stoichiometric_dict["CDGUNPD"]
        self.stoichiometric_dict["CDPDSP_BS"] = {}
        self.stoichiometric_dict["G3POA_BS"] = {}
        self.stoichiometric_dict["G3PCT"] = self.stoichiometric_dict["G3PCtex"]
        self.stoichiometric_dict["LYSLG_BS"] = {}
        self.stoichiometric_dict["KAS1"] = self.stoichiometric_dict["KAS13"]
        self.stoichiometric_dict["PGPPH_BS"] = {}
        self.stoichiometric_dict["PPTGS_BS"] = {}
        self.stoichiometric_dict["PSDC_BS"] = {}
        self.stoichiometric_dict["TECA1S45"] = {}
        self.stoichiometric_dict["TECA4S_BS"] = {}
        self.stoichiometric_dict["UAG2E"] = self.stoichiometric_dict["UAG2EMA"]
        self.stoichiometric_dict["UAG4E"] = self.stoichiometric_dict["UAG4Ei"]
        self.stoichiometric_dict["UGT1_BS"] = {}
        self.stoichiometric_dict["UGT2_BS"] = {}
        self.stoichiometric_dict["UGT_BS"] = {}

        self.active_network_stoich = {key: self.stoichiometric_dict[key]
                                      for key in self.list_of_active_reactions}
        for reaction_id, reaction_stoich in self.active_network_stoich.items():
            self.active_reactions.append(reaction_id)
            for metabolite_id, stoich_coeff in reaction_stoich.items():
                if metabolite_id not in self.active_metabolites:
                    self.active_metabolites.append(metabolite_id)
            # print(reaction_id, reaction_stoich)

    def build_cobra_model(self):
        """
        This method, build a metabolic network in the COBRA format based on self.active_metabolites,
            self.active_reactions, and self.active_network_stoich.
        :return: -, filling self.cobra_model.
        """
        model = Model('active_model')
        model_metabolites = [Metabolite(met_id) for met_id in self.active_metabolites]
        model.add_metabolites(model_metabolites)
        model_reactions = [Reaction(rxn_id) for rxn_id in self.active_reactions]
        model.add_reactions(model_reactions)
        for model_reaction in model.reactions:
            model_reaction.add_metabolites(self.active_network_stoich[model_reaction.id])
        self.cobra_model = model

    def visualize_active_network(self, filepath_to_save_html: str):
        """
        This method, makes the plot for the self.cobra_model, and saves it in filepath_to_save_html.
        :param filepath_to_save_html: Filepath to save .html plot
        :return: -
        """
        # graph = nx.DiGraph()
        # pos = {}
        # edgelist = []
        # for met_index in range(len(self.active_metabolites)):
        #     met_id = self.active_metabolites[met_index]
        #     pos[met_id] = (0.25, met_index * 0.1)
        # for rxn_index in range(len(self.active_reactions)):
        #     rxn_id = self.active_reactions[rxn_index]
        #     pos[rxn_id] = (0.75, rxn_index * 0.1)
        # nx.draw_networkx_nodes(graph, pos, node_size=10, nodelist=self.active_metabolites, node_color="tab:blue")
        # nx.draw_networkx_nodes(graph, pos, node_size=10, nodelist=self.active_reactions, node_color="tab:red")
        # for rxn_id, rxn_mets in self.active_network_stoich.items():
        #     for met_id, met_coeff in rxn_mets.items():
        #         if met_coeff > 0:
        #             edgelist.append((rxn_id, met_id))
        #         else:
        #             edgelist.append((met_id, rxn_id))
        # nx.draw_networkx_edges(graph, pos, edgelist=edgelist, edge_color="tab:green")
        # plt.show()
        self.build_cobra_model()
        # ecoli_model = cobra.io.read_sbml_model("../Data/Palsson B.Subtilis Reconstruction/Results/e_coli.xml")
        builder = escher.Builder(map_json="../Data/Palsson B.Subtilis Reconstruction/Results/central_metabolism.json",
                                 model=self.cobra_model)
        # builder = escher.Builder(map_json="../Data/Palsson B.Subtilis Reconstruction/Results/central_metabolism.json",
        #                          model=ecoli_model)
        builder.save_html(filepath_to_save_html)


final_V = pd.read_csv("../Data/Palsson B.Subtilis Reconstruction/Results/final_csv.csv")
reactions_fluxes = final_V['x210']
# reactions_activation = reaction_fluxes.any(axis='columns')
active_reactions_indexes = reactions_fluxes.to_numpy().nonzero()[0]

rxn_map_path = "../Data/Palsson B.Subtilis Reconstruction/Microbial Final Data/reactions_index_map.json"
with open(rxn_map_path, 'r') as js_file:
    rxn_map = json.load(js_file)
rxn_rev_map = {value: key for (key, value) in rxn_map.items()}
active_reactions_ids = [rxn_rev_map[rxn_id] for rxn_id in active_reactions_indexes]

compare_rpi_df = pd.read_csv("../Data/Palsson B.Subtilis Reconstruction/Results/rpi_nz.csv")
compare_rpi_reactions = compare_rpi_df[compare_rpi_df.columns[0]].tolist()

stoich_filepath = "../Data/Palsson B.Subtilis Reconstruction/Microbial Template/Microbial Stoichiometric Data.json"
net_visualizer = ActiveNetworkVisualizer(list_of_active_reactions=active_reactions_ids,
                                         stoichiometry_dict_filepath=stoich_filepath)
html_saving_filepath = "../Data/Palsson B.Subtilis Reconstruction/Results/result_net.html"
net_visualizer.visualize_active_network(filepath_to_save_html=html_saving_filepath)

rpi_visualizer = ActiveNetworkVisualizer(list_of_active_reactions=compare_rpi_reactions,
                                         stoichiometry_dict_filepath=stoich_filepath)
html_saving_filepath = "../Data/Palsson B.Subtilis Reconstruction/Results/rpi_net.html"
rpi_visualizer.visualize_active_network(filepath_to_save_html=html_saving_filepath)
