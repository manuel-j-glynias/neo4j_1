import sys

import requests
import json
import csv
import uniprot
from neo4j import GraphDatabase

from NCI_Drugs import get_nci_target_list


def get_io_list():
    # file_path = 'data/tblPROD_OA_OA_MarkerReporting.csv'
    # file_path = 'data/1_MarkerReporting.csv'
    file_path = 'data/Small_MarkerReporting.csv'
    io_list = []
    pheno_list = []
    with open(file_path, newline='') as io_data:
        io_reader = csv.DictReader(io_data)
        for io in io_reader:
            if len(io['Therapies']) > 0:
                io['drug_list'] = io['Therapies'].split(";")
            if io['ImmunePhenotype'] not in pheno_list:
                pheno_list.append(io['ImmunePhenotype'])

            io_list.append(dict(io))

    return io_list, pheno_list



def get_io_and_pheno_lists():
    io_list, pheno_list = get_io_list()
    for io in io_list:
        gene_info, sp_info = uniprot.add_gene_and_sp_info(io['GeneSymbol'])
        io['gene_info'] = gene_info
        io['sp_info'] = sp_info
        # print(io['GeneSymbol'], io['sp_info']['transcript'])
    return io_list, pheno_list


def send_mutation(mutation_payload):
    url = "http://localhost:7474/graphql/"
    headers = {
      'Authorization': 'Basic bmVvNGo6b21uaQ==',
      'Content-Type': 'application/json',
    }

    response = requests.request("POST", url, headers=headers, data = mutation_payload)
    if not response.ok:
        response.raise_for_status()
        sys.exit()

    responseBody = response.json()

    return responseBody


def get_mutation_start(muatation_number):
    mutation_id = f'm{muatation_number}'
    s = 'mutation ' + str(mutation_id) + ' {\n'
    return s

def write_phenotype_mutation(pheno_list, pheno_dict):
    mutation_payload = '{"query":"mutation {'
    s = ''
    for index, pheno in enumerate(pheno_list, start=1):
        id = 'ph' + str(index)
        pheno_dict[pheno] = id
        if len(s) > 0:
            s += ','
        s += f'{id}: createImmunePhenotype(id: \\"{id}\\", name: \\"{pheno}\\")'
    mutation_payload += s
    mutation_payload += '}"}'
    responseBody = send_mutation(mutation_payload)
    print(responseBody)


def get_name_for_internet_reference(url,accessed_date):
    pos1 = url.find('//') + 2
    pos2 = url.find('/',pos1)
    name = url[pos1:pos2]
    name += ' (accessed on:' + accessed_date + ')'
    return name


def get_acessed_date_as_string(d):
    date_time = d.strftime("%m/%d/%Y")
    return date_time


def write_gene_mutation(index, component):
    id = 'g' + str(index)
    gene_info = component['gene_info']
    gene = component['GeneSymbol']
    chrom = gene_info['chrom']
    strand = gene_info['strand']
    start = gene_info['start']
    end = gene_info['end']
    statement = gene + ' on chromosome ' + chrom + ' at ' + str(start) + '-' + str(end)
    s = f'{id}: createGene(id: \\"{id}\\", name: \\"{gene}\\", chromosome: \\"{chrom}\\", strand:{strand}, start: {start}, end: {end}, statement:\\"{statement}\\"),'
    ref = gene_info['reference']
    if ref['type']=='InternetReference':
        ref_id = 'ref_' + id
        accessed = get_acessed_date_as_string(ref['accessed_date'])
        ir_name = get_name_for_internet_reference(ref["url"],accessed)
        s += f'{ref_id}: createInternetReference(id:\\"{ref_id}\\", shortReference: \\"{ir_name}\\", web_address:\\"{ref["url"]}\\", accessed_date:\\"{accessed}\\"),'
        ref_id2 = 'gref_' + id
        s += f'{ref_id2}: addGeneReference(id:\\"{id}\\", reference:\\"{ref_id}\\" ),'
    return id,s

def write_polypeptide_mutation(index, component, gene_id, isoform):
    sp_info = component['sp_info']
    id = 'pp' + str(index)
    polypeptide = sp_info['name']
    uniprot_id = sp_info['id']
    acc_num = sp_info['acc_num']
    s = f'{id}: createPolypeptide(id: \\"{id}\\", name: \\"{polypeptide}\\", uniprot_id: \\"{uniprot_id}\\", accessionNumber: \\"{acc_num}\\" ),'

    # AddPolypeptideGene
    tag = id + '_' + gene_id
    s += f'{tag}: addPolypeptideGene(gene:\\"{gene_id}\\", id:\\"{id}\\"),'

    # AddPolypeptideTranscript
    tag = id + '_' + isoform
    s += f'{tag}: addPolypeptideTranscript(transcript:\\"{isoform}\\", id:\\"{id}\\"),'
    return id,s


def write_isoform_mutation(index, component,gene_id):
    id = 'iso' + str(index)
    isoform = component['sp_info']['transcript']
    s = f'{id}: createIsoform(id: \\"{id}\\", name: \\"{isoform}\\", ),'

    tag = id + '_' + gene_id

    s += f'{tag}: addIsoformGene(gene:\\"{gene_id}\\",id:\\"{id}\\"),'
    return id, s


def write_protein_mutation(index, io):
    s = ''
    id = 'protein_' + str(index)
    protein = io['ProteinMarker']
    s = f'{id}: createProtein(id: \\"{id}\\", name: \\"{protein}\\", synonyms: []),'
    return id,s



def write_protein_polypeptide_mutation(protein, polypeptides):
    s = ''
    tag = protein
    pp_array = '['
    for polypeptide in polypeptides:
        tag += '_' + polypeptide
        pp_array += f'\\"{polypeptide}\\"'
    pp_array += ']'
    s = f'{tag}: addProteinPolypeptide_chains(id:\\"{protein}\\", polypeptide_chains:{pp_array}),'
    return s


def write_protein_phenotype_mutation(protein, phenotype):
    tag = protein + '_' + phenotype
    s = f'{tag}: addImmunePhenotypeBiomarkers(id:\\"{phenotype}\\", proteins:[\\"{protein}\\"]),'
    return s


def add_createSubcellularLocationWithEvidence(sclwe_id,statement):
    s = f'{sclwe_id}: createSubcellularLocationWithEvidence(id:\\"{sclwe_id}\\", statement:\\"{statement}\\"),'
    return s

def create_addsubcellularLocationWithEvidenceLocation_mutation(sclwe_id, location_id):
    id = sclwe_id + '_' + location_id
    s = f'{id}: addSubcellularLocationWithEvidenceLocation(id:\\"{sclwe_id}\\", location:[\\"{location_id}\\"]),'
    return s

def create_addsubcellularLocationWithEvidenceReference_mutation(sclwe_id, ref_ids):
    id = 'ref_' + sclwe_id
    reference_string = '['
    for r in ref_ids:
        if len(reference_string)>1:
            reference_string += ","
        reference_string += '\\"' + r + '\\"'
    reference_string += ']'

    s = f'{id}: addSubcellularLocationWithEvidenceReference(id:\\"{sclwe_id}\\", reference:{reference_string}),'
    return s


def create_addProteinLocation(protein_id, sclwe_id):
    id = sclwe_id + '_' + protein_id
    s = f'{id}: addProteinLocation(id:\\"{protein_id}\\", location:\\"{sclwe_id}\\"),'
    return s


def add_createSubcellularLocation_mutation(location_id, location):
    s = f'{location_id}: createSubcellularLocation(id:\\"{location_id}\\", name:\\"{location}\\"),'
    return s


def create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id):
    id = ref_id + '_' + journal_id
    s = f'{id}: addLiteratureReferenceJournal(id:\\"{ref_id}\\", journal:\\"{journal_id}\\"),'
    return s



def create_AddLiteratureReferenceAuthors_mutation(ref_id, authors):
    id = 'au_' +ref_id
    author_string = '['
    for a in authors:
        if len(author_string)>1:
            author_string += ","
        author_string += '\\"' + a + '\\"'
    author_string += ']'
    s = f'{id}: addLiteratureReferenceAuthors(id:\\"{ref_id}\\", authors:{author_string}),'

    return s


def create_reference_mutation(ref_id, ref):
    ref_name = ref_name_from_authors_pmid_and_year(ref['authors'], ref['PubMed'], ref['year'])
    s = f'''{ref_id}: createLiteratureReference(id: \\"{ref_id}\\", shortReference: \\"{ref_name}\\", title: \\"{ref['title']}\\", volume: \\"{ref['volume']}\\", first_page: \\"{ref['first_page']}\\", last_page: \\"{ref['last_page']}\\", publication_Year: \\"{ref['year']}\\", DOI: \\"{ref['DOI']}\\", PMID: \\"{ref['PubMed']},\\"),'''
    return s


def create_author_mutation(id,surname,first,middle):
    s = f'''{id}: createAuthor(id: \\"{id}\\",surname: \\"{surname}\\",first_initial: \\"{first}\\" ,middle_initial: \\"{middle}\\"),'''
    return s


def create_journal_mutation(journal, journal_id):
    s = f'''{journal_id}: createJournal(id: \\"{journal_id}\\",name: \\"{journal}\\"),'''
    return s


def ref_name_from_authors_pmid_and_year(authors, pmid, year):
    s = ''
    first, middle, surname = get_authors_names(authors[0])
    if len(authors) == 1:
        s += surname + ' ' + year
    elif len(authors) == 2:
        first2, middle2, surname2 = get_authors_names(authors[1])
        s += surname + '& '+ surname2 + ' ' + year
    else:
        s += surname + ' et al. ' + year
    s += ' (PMID:' + pmid + ')'
    return s


def write_location_mutations(polypeptide_info,author_dict,journal_dict,reference_dict,location_dict,protein_id,protein_name):
    s = ''
    location = polypeptide_info['location']
    references = polypeptide_info['location_references']
    journal_id = None
    ref_id = None

    if location in location_dict:
        location_id = location_dict[location]
    else:
        location_id = 'loc_' + str(len(location_dict) + 1)
        location_dict[location] = location_id
        s += add_createSubcellularLocation_mutation(location_id,location)
    refs = []
    for ref in references:
        if ref['type'] == 'LiteratureReference':
            pubmed = ref['PubMed']
            if not pubmed in reference_dict:
                ref_id = 'ref_' + str(len(reference_dict) + 1)
                reference_dict[pubmed] = ref_id
                s += create_reference_mutation(ref_id,ref)
                journal = ref['journal']
                if not journal in journal_dict:
                    journal_id = 'j_' + str(len(journal_dict) + 1)
                    journal_dict[journal] = journal_id
                    s += create_journal_mutation(journal, journal_id)
                else:
                    journal_id = journal_dict[journal]
                s += create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id)
                authors = []
                for author in ref['authors']:
                    first, middle, surname = get_authors_names(author)
                    id = 'au_' + surname + '_' + first + '_' + middle
                    id = id.replace("-", "_")
                    if not id in author_dict:
                        author_dict[id] = id
                        s += create_author_mutation(id, surname, first, middle)
                        authors.append(id)
                s += create_AddLiteratureReferenceAuthors_mutation(ref_id, authors)

            else:
                ref_id = reference_dict[pubmed]
            refs.append(ref_id)
    sclwe_id = 'sclwe_' + location_id + '_' + protein_id
    statement = protein_name + ' at ' + location
    s += add_createSubcellularLocationWithEvidence(sclwe_id,statement)
    s += create_addsubcellularLocationWithEvidenceLocation_mutation(sclwe_id, location_id)
    s += create_addsubcellularLocationWithEvidenceReference_mutation(sclwe_id, refs)
    s += create_addProteinLocation(protein_id, sclwe_id)
    return s


def get_authors_names(author):
    l = author.split()
    surname = l[0].replace("'", "_")
    ll = l[1].split('.')
    first = ll[0]
    if len(ll) > 1 and ll[1] is not '':
        middle = ll[1]
    else:
        middle = '_'
    return first, middle, surname

def return_graphql_boolean(b):
    s = 'false'
    if b:
        s = 'true'
    return s


def write_variant_drugs(drug_list, drug_dict, variant_id):
    s = ''
    if len(drug_list) > 0:
        drug_id_string = '['
        for drug in drug_list:
            drug_name = drug['PreferredName'].strip()
            if drug_name in drug_dict:
                drug_id = drug_dict[drug_name]
            else:
                drug_id = 'd_' + str(len(drug_dict) + 1)
                drug_dict[drug_name] = drug_id
                synonyms = '['
                for syn in drug['synonyms']:
                    synonyms += f'\\"{syn}\\",'
                synonyms += ']'
                s += f'''{drug_id}: createDrug(annotationDate:\\"{drug['AnnotationDate']}\\",concept:\\"{drug['ConceptCode']}\\",definition:\\"{drug['Definition']}\\",drugCategory:\\"{drug['DrugCategory']}\\",id: \\"{drug_id}\\", isAntineoplastic:{return_graphql_boolean(drug['Antineoplastic_Flag'])},isImmuno:{return_graphql_boolean(drug['Immuno_Flag'])}, modulator:\\"{drug['Modulator']}\\", synonyms:{synonyms}, name: \\"{drug_name}\\"),'''
            if len(drug_id_string) > 1:
                drug_id_string += ","
            drug_id_string += '\\"' + drug_id + '\\"'
        drug_id_string += ']'
        id = variant_id + '_drugs'
        s += f'{id}: addVariantTargetDrugs(id:\\"{variant_id}\\", drugs:{drug_id_string}),'
    return s


def write_protein_target_drugs(drug_list, drug_dict, protein_target_id):
    s = ''
    if len(drug_list) > 0:
        drug_id_string = '['
        for drug in drug_list:
            drug_name = drug['PreferredName'].strip()
            if drug_name in drug_dict:
                drug_id = drug_dict[drug_name]
            else:
                drug_id = 'd_' + str(len(drug_dict) + 1)
                drug_dict[drug_name] = drug_id
                synonyms = '['
                for syn in drug['synonyms']:
                    synonyms += f'\\"{syn}\\",'
                synonyms += ']'
                s += f'''{drug_id}: createDrug(annotationDate:\\"{drug['AnnotationDate']}\\",concept:\\"{drug['ConceptCode']}\\",definition:\\"{drug['Definition']}\\",drugCategory:\\"{drug['DrugCategory']}\\",id: \\"{drug_id}\\", isAntineoplastic:{return_graphql_boolean(drug['Antineoplastic_Flag'])},isImmuno:{return_graphql_boolean(drug['Immuno_Flag'])}, modulator:\\"{drug['Modulator']}\\", synonyms:{synonyms}, name: \\"{drug_name}\\"),'''
                # s += f'''{drug_id}: createDrug(annotationDate:\\"{drug['AnnotationDate']}\\",concept:\\"{drug['ConceptCode']}\\",definition:\\"{drug['Definition']}\\",drugCategory:\\"{drug['DrugCategory']}\\",id: \\"{drug_id}\\", isAntineoplastic:{return_graphql_boolean(drug['Antineoplastic_Flag'])},isImmuno:{return_graphql_boolean(drug['Immuno_Flag'])}, modulator:\\"{drug['Modulator']}\\", , synonyms:[], name: \\"{drug_name}\\"),'''

            if len(drug_id_string) > 1:
                drug_id_string += ","
            drug_id_string += '\\"' + drug_id + '\\"'

        drug_id_string += ']'
        id = protein_target_id + '_drugs'
        s += f'{id}: addProteinTargetDrugs(id:\\"{protein_target_id}\\", drugs:{drug_id_string}),'
    return s

def write_ProteinTarget_mutation(marker_object, protein,drug_dict):
    protein_target_id = 'protTarget_' + marker_object['ProteinMarker']
    name = marker_object['ProteinMarker'] + ' wildtype'
    s = f'''{protein_target_id}: createProteinTarget(id:\\"{protein_target_id}\\",name:\\"{name}\\"),'''
    s += f'''addProteinTargetProtein(id:\\"{protein_target_id}\\", protein:[\\"{protein}\\"]),'''
    if 'drug_list' in marker_object:
        m = write_protein_target_drugs(marker_object['drug_list'], drug_dict, protein_target_id)
        s += m
    return s

def write_VariantTarget_mutation(marker_object, protein,drug_dict,var_counter):
    variant_id = 'varTarget_' + str(var_counter)
    s = f'''{variant_id}: createVariantTarget(id:\\"{variant_id}\\",name:\\"{marker_object['variantTarget']}\\"),'''
    mutation_id = variant_id+ '_' + protein
    s += f'''{mutation_id}: addVariantTargetProtein(id:\\"{variant_id}\\", protein:[\\"{protein}\\"]),'''
    if 'drug_list' in marker_object:
        m = write_variant_drugs(marker_object['drug_list'], drug_dict, variant_id)
        s += m
    return s

def write_pathways_mutations(pathways, pathways_dict, protein):
    s = ''
    pathway_ids = []
    for pathway in pathways:
        if pathway in pathways_dict:
            pathway_id = pathways_dict[pathway]
        else:
            pathway_id = 'pw_' + str(len(pathways_dict) + 1)
            pathways_dict[pathway] = pathway_id
            s += f'''{pathway_id}: createPathway(id: \\"{pathway_id}\\",name: \\"{pathway}\\"),'''
        pathway_ids.append(pathway_id)
    if len(pathway_ids)> 0:
        pathway_id_string = '['
        for pw in pathway_ids:
            if len(pathway_id_string) > 1:
                pathway_id_string += ","
            pathway_id_string += '\\"' + pw + '\\"'
        pathway_id_string += ']'
        id = protein + '_pathways'
        s += f'{id}: addProteinPathways(id:\\"{protein}\\", pathways:{pathway_id_string}),'
    return s




def write_mutation(marker_list, pheno_list):
    pheno_dict = {}
    protein_dict = {}
    gene_dict = {}
    polypeptide_dict = {}
    author_dict = {}
    journal_dict = {}
    reference_dict = {}
    location_dict = {}
    drug_dict = {}
    pathways_dict = {}

    gene_counter = 0
    var_counter = 0
    if pheno_list is not None:
        write_phenotype_mutation(pheno_list,pheno_dict)
    for index, marker_object in enumerate(marker_list, start=1):
        pp_array = []
        mutation_payload = '{"query":"mutation {'
        s = ''
        if marker_object['ProteinMarker'] in protein_dict:
            protein = protein_dict[marker_object['ProteinMarker']]
        else:
            for component in marker_object['components']:
                if component['GeneSymbol'] in gene_dict:
                    gene_id = gene_dict[component['GeneSymbol']]
                else:
                    gene_id,m = write_gene_mutation(gene_counter, component)
                    s += m
                    gene_dict[component['GeneSymbol']] = gene_id
                    isoform,m = write_isoform_mutation(gene_counter, component,gene_id)
                    s += m
                polypeptide_name = component['sp_info']['name']
                if polypeptide_name in polypeptide_dict:
                    polypeptide = polypeptide_dict[polypeptide_name]
                else:
                    polypeptide,m = write_polypeptide_mutation(gene_counter, component,gene_id,isoform)
                    s += m
                    polypeptide_dict[polypeptide_name] = polypeptide
                pp_array.append(polypeptide)
                gene_counter += 1
            protein,m = write_protein_mutation(index,marker_object)
            s += m
            protein_dict[marker_object['ProteinMarker']] = protein
            m = write_protein_polypeptide_mutation(protein,pp_array)
            s += m
            m = write_location_mutations(marker_object['components'][0]['sp_info'], author_dict, journal_dict, reference_dict, location_dict,
                                     protein, marker_object['ProteinMarker'])
            s += m
            m = write_pathways_mutations(marker_object['components'][0]['sp_info']['pathways'], pathways_dict, protein)
            s += m
            if pheno_list is not None:
                phenotype = pheno_dict[marker_object['ImmunePhenotype']]
                m = write_protein_phenotype_mutation(protein, phenotype)
                s += m
        if marker_object['variantTarget']==None:
            m = write_ProteinTarget_mutation(marker_object,protein,drug_dict)
            s += m
        else:
            m = write_VariantTarget_mutation(marker_object,protein,drug_dict,var_counter)
            s += m
            var_counter += 1
        s += '}'
        mutation_payload += s
        mutation_payload += '"}'
        print(mutation_payload)
        responseBody = send_mutation(mutation_payload)
        print(responseBody)



def main():
    uri = "bolt://localhost:7687"

    with open('schema.graphql', 'r') as file:
        idl_as_string = file.read()

    driver = GraphDatabase.driver(uri, auth=("neo4j", "omni"))
    with driver.session() as session:
        tx = session.begin_transaction()
        tx.run("match(a) detach delete(a)")
        result = tx.run("call graphql.idl('" + idl_as_string + "')")
        print(result.single()[0])
        tx.commit()
    driver.close()

    g_marker_list = ['CD8','EGFR', 'BRAF', 'MET']
    # g_marker_list = ['EGFR', 'BRAF', 'MET']
    gene_dict = {}
    protein_dict = {}
    targets = get_nci_target_list(g_marker_list)
    for target in targets:
        uniprot.add_gene_and_sp_info(target, gene_dict, protein_dict)
    write_mutation(targets, None)



if __name__ == "__main__":
    main()

