from GenesKO_Standardizer import GenesKOStandardizer
from GPR_MapStandardizer import GPRMapConverter
from ReactionsKOMaker import ReactionsKOMaker
from KnockOutBoundsMaker import KnockOutBoundsMaker
from SourceUtilStandardizer import SourceUtilGrowthData
from SourceUtilBoundsMaker import SourceUtilBoundsMaker
from TemplateBoundsMaker import TemplateBoundsMaker


class BioDataOrganizer:  # ToDo: 1.descriptions  2.run from middle (we have mid data)
    def __init__(self,
                 media_filepaths_dict: dict,
                 internal_rxns_filepath: str,
                 folder_to_save_final_bounds: str = "./"):
        self.media_filepaths_dict = media_filepaths_dict
        self.internal_rxns_filepath = internal_rxns_filepath
        self.folder_to_save_final_bounds = folder_to_save_final_bounds

    def organize_ko_bounds(self,
                           genes_ko_growth_filepath: str,
                           medium_name: str,
                           gene_assoc_data_filepath: str,
                           input_reactions_nomenclature: str,
                           output_reactions_nomenclature: str,
                           translation_filepath: str,
                           filepath_to_save_ko_genes_dict: str = "./Genes KO Growth.json",
                           filepath_to_save_organism_gpr: str = "./Organism GPR.json",
                           filepath_to_save_ko_reactions_dict: str = "./Reactions KO Growth.json"):
        # ################### Convert KO data into the standard .json ####################
        genes_ko_standardizer = GenesKOStandardizer(
            genes_ko_growth_filepath=genes_ko_growth_filepath,
            medium_name=medium_name)
        genes_ko_standardizer.make_genes_ko_growth_dict(
            filepath_to_save=filepath_to_save_ko_genes_dict)
        # ################# Convert GPR association data into the standard .json ##################
        gpr_map_converter = GPRMapConverter(
            gene_assoc_data_filepath=gene_assoc_data_filepath,
            input_reactions_nomenclature=input_reactions_nomenclature,
            output_reactions_nomenclature=output_reactions_nomenclature,
            translation_filepath=translation_filepath)
        gpr_map_converter.make_genes_to_reactions_ko_dict(
            filepath_to_save=filepath_to_save_organism_gpr)
        # ##################### Convert genes KO data into reactions KO data ######################
        reaction_ko_maker = ReactionsKOMaker(
            organism_gpr_filepath=filepath_to_save_organism_gpr,
            genes_ko_growth_filepath=filepath_to_save_ko_genes_dict)
        reaction_ko_maker.make_reactions_ko_growth(
            filepath_to_save=filepath_to_save_ko_reactions_dict)
        # ######################### Making and saving final KO bounds ##########################
        ko_bounds_maker = KnockOutBoundsMaker(
            reactions_ko_filepath=filepath_to_save_ko_reactions_dict,
            media_filepath_dict=self.media_filepaths_dict,
            internal_rxns_filepath=self.internal_rxns_filepath)
        ko_bounds_maker.make_growth_bounds()
        ko_bounds_maker.make_non_growth_bounds()
        ko_bounds_maker.save_all_bounds(folder_to_save=self.folder_to_save_final_bounds + 'KO Bounds/')

    def organize_source_util_bounds(self,
                                    source_util_csv_filepath: str,
                                    medium_name: str,
                                    filepath_to_save_util_dict: str = "./Growth_Biolog.json"):
        # ################### Convert util. data into the standard .json ####################
        source_util_data = SourceUtilGrowthData(
            csv_filepath=source_util_csv_filepath,
            medium_name=medium_name)
        source_util_data.make_json_file(
            filepath_to_save=filepath_to_save_util_dict)
        # ################### Making and saving final source util. bounds ####################
        util_bounds_maker = SourceUtilBoundsMaker(
            sources_util_filepath=filepath_to_save_util_dict,
            media_filepath_dict=self.media_filepaths_dict,
            internal_rxns_filepath=self.internal_rxns_filepath)
        util_bounds_maker.make_growth_bounds()
        util_bounds_maker.make_non_growth_bounds()
        util_bounds_maker.save_all_bounds(folder_to_save=self.folder_to_save_final_bounds + 'Util Bounds/')

    def finalize_bounds(self):
        template_bounds_maker = TemplateBoundsMaker([], [])
        print(template_bounds_maker)


media_dict = {'minimal_media': "../Data/Palsson B.Subtilis Reconstruction/Biolog_Medium_Bounds.csv",
              'LB_Rich_Medium': "../Data/Palsson B.Subtilis Reconstruction/LB_Medium_Bounds.csv"}
bio_data_organizer = BioDataOrganizer(
    media_filepaths_dict=media_dict,
    internal_rxns_filepath="../Data/Palsson B.Subtilis Reconstruction/Internal_Rxns_Bounds.csv",
    folder_to_save_final_bounds="../Data/Palsson B.Subtilis Reconstruction/")

bio_data_organizer.organize_ko_bounds(
    genes_ko_growth_filepath="../Data/Palsson B.Subtilis Reconstruction/Genes KO Growth Data.csv",
    medium_name="LB_Rich_Medium",
    filepath_to_save_ko_genes_dict="../Data/Palsson B.Subtilis Reconstruction/Genes KO Growth.json",
    gene_assoc_data_filepath="../Data/Palsson B.Subtilis Reconstruction/Genes Associations.json",
    input_reactions_nomenclature='R_Rxn',
    output_reactions_nomenclature='Rxn',
    translation_filepath="../Data/Palsson B.Subtilis Reconstruction/Reactions Translation File.csv",
    filepath_to_save_organism_gpr="../Data/Palsson B.Subtilis Reconstruction/Organism GPR.json",
    filepath_to_save_ko_reactions_dict="../Data/Palsson B.Subtilis Reconstruction/Reactions KO Growth.json")

bio_data_organizer.organize_source_util_bounds(
    source_util_csv_filepath="../Data/Palsson B.Subtilis Reconstruction/Growth_Biolog.csv",
    medium_name="minimal_media",
    filepath_to_save_util_dict="../Data/Palsson B.Subtilis Reconstruction/Growth_Biolog.json")
