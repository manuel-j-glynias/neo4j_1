import csv
import pprint


def get_concept_code_dict():
    file_path = 'data/tblOS_GLOBAL_NCI_DL_Concepts.csv'
    drug_dict = {}
    with open(file_path, newline='', encoding='utf-8',errors='ignore') as drug_data:
        drug_reader = csv.DictReader(drug_data)
        for drug in drug_reader:
            drug_dict[drug['ConceptCode']] = drug
    return drug_dict


def get_synonym_dict():
    file_path = 'data/tblOS_GLOBAL_NCI_DL_Synonyms.csv'
    synonym_dict = {}
    with open(file_path, newline='', encoding='utf-8',errors='ignore') as drug_data:
        drug_reader = csv.DictReader(drug_data)
        for drug in drug_reader:

            synonym = drug['Synonym'].replace("'","`")
            if drug['ConceptCode'] in synonym_dict:
                synonym_dict[drug['ConceptCode']] ['synonyms'].append(synonym)
            else:
                concept = {'ConceptCode':drug['ConceptCode'], 'synonyms':[synonym]}
                synonym_dict[drug['ConceptCode']] = concept
    return synonym_dict


def get_subset_of_nci_drugs_and_targets(g_marker_list,drug_dict,synonym_dict):
    file_path = 'data/tblOS_GLOBAL_NCI_DL_Targets.csv'
    target_list = []
    target_dict = {}
    with open(file_path, newline='') as drug_data:
        drug_reader = csv.DictReader(drug_data)
        for drug in drug_reader:
            if drug['Skip_Flag'] == '0':
                g_marker = drug['HGNC_Symbol'].split(' ')[0]
                if g_marker_list is None or g_marker in g_marker_list:
                    concept = drug_dict[drug['ConceptCode']]
                    synonyms = []
                    if drug['ConceptCode'] in synonym_dict:
                        synonyms = synonym_dict[drug['ConceptCode']]['synonyms']
                    drug_object = {'ConceptCode': drug['ConceptCode'], 'PreferredName': drug['PreferredName'],
                                   'AnnotationDate': drug['AnnotationDate'],'synonyms':synonyms,
                                   'Antineoplastic_Flag': concept['Antineoplastic_Flag']=='-1','Immuno_Flag':concept['Immuno_Flag']=='-1',
                                   'DrugCategory':concept['DrugCategory'],'Definition':concept['Definition'],'Modulator':concept['Modulator']}

                    variant_string = drug['GeneVariant']
                    variant_string = variant_string.replace(' and/or ',',')
                    variant_string = variant_string.replace(' and ',',')
                    variants = variant_string.split(',')
                    for variant in variants:
                        variant = variant.strip()
                        if len(variant)>0 and not variant.startswith(g_marker):
                            variant = g_marker + ' ' + variant
                        variantTarget = None
                        key = g_marker
                        if len(variant)>0 and 'wildtype' not in variant:
                            variantTarget = variant
                            key = variant

                        if key in target_dict:
                            target_dict[key]['drug_list'].append(drug_object)
                        else:
                            gene_object = {'ProteinMarker':g_marker,'GeneSymbol': g_marker, 'drug_list':[drug_object], 'variantTarget':variantTarget}
                            target_list.append(gene_object)
                            target_dict[key] = gene_object
    return target_list


def main():
    target_list = get_nci_target_list()
    pprint.pprint(target_list)


def get_nci_target_list(g_marker_list):
    drug_dict = get_concept_code_dict()
    synonym_dict = get_synonym_dict()
    target_list = get_subset_of_nci_drugs_and_targets(g_marker_list, drug_dict, synonym_dict)
    return target_list


if __name__ == "__main__":
    main()

