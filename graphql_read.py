
import json

import csv
import uniprot

def get_io_list():
    file_path = 'data/tblPROD_OA_OA_MarkerReporting.csv'
    # file_path = 'data/1_MarkerReporting.csv'
    # file_path = 'data/Small_MarkerReporting.csv'
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



def write_location_mutations(protein_info,author_dict,journal_dict,reference_dict,location_dict,biomarker):
    s = ''
    location = protein_info['location']
    references = protein_info['location_references']
    journal_id = None
    ref_id = None

    if location in location_dict:
        location_id = location_dict[location]
    else:
        location_id = 'loc_' + str(len(location_dict) + 1)
        location_dict[location] = location_id
        s += add_createSubcellularLocation_mutation(location_id,location)

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
                    l = author.split()
                    surname = l[0].replace("'", "_")

                    ll = l[1].split('.')
                    first = ll[0]
                    if len(ll) > 1 and ll[1] is not '':
                        middle = ll[1]
                    else:
                        middle = '_'
                    id = 'au_' + surname + '_' + first + '_' + middle
                    id = id.replace("-", "_")
                    if not id in author_dict:
                        author_dict[id] = id
                        s += create_author_mutation(id, surname, first, middle)
                        authors.append(id)
                s += create_AddLiteratureReferenceAuthors_mutation(ref_id, authors)

            else:
                ref_id = reference_dict[pubmed]
            sclwe_id = location_id + '_' + ref_id
            s += add_createSubcellularLocationWithEvidence(sclwe_id)
            s += create_addsubcellularLocationWithEvidenceLocation_mutation(sclwe_id, location_id)
            s += create_addsubcellularLocationWithEvidenceReference_mutation(sclwe_id, ref_id)
            s += create_addBiomarkerLocation(biomarker,sclwe_id)
    return s

def add_createSubcellularLocationWithEvidence(sclwe_id):
    s = f'{sclwe_id}: createSubcellularLocationWithEvidence(id:"{sclwe_id}")\n'
    return s

def create_addsubcellularLocationWithEvidenceLocation_mutation(sclwe_id, location_id):
    id = sclwe_id + '_' + location_id
    s = f'{id}: addSubcellularLocationWithEvidenceLocation(id:"{sclwe_id}", location:"{location_id}")\n'
    return s

def create_addsubcellularLocationWithEvidenceReference_mutation(sclwe_id, ref_id):
    id = sclwe_id + '_' + ref_id
    s = f'{id}: addSubcellularLocationWithEvidenceReference(id:"{sclwe_id}", reference:"{ref_id}")\n'
    return s


def create_addBiomarkerLocation(biomarker_id,sclwe_id):
    id = sclwe_id + '_' + biomarker_id
    s = f'{id}: addBioMarkerLocation(id:"{biomarker_id}", location:"{sclwe_id}")\n'
    return s


def add_createSubcellularLocation_mutation(location_id, location):
    s = f'{location_id}: createSubcellularLocation(id:"{location_id}", name:"{location}")\n'
    return s


def create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id):
    id = ref_id + '_' + journal_id
    s = f'{id}: addLiteratureReferenceJournal(id:"{ref_id}", journal:"{journal_id}")\n'
    return s



def create_AddLiteratureReferenceAuthors_mutation(ref_id, authors):
    id = 'au_' +ref_id
    author_string = '['
    for a in authors:
        if len(author_string)>1:
            author_string += ","
        author_string += '"' + a + '"'
    author_string += ']'
    s = f'{id}: addLiteratureReferenceAuthors(id:"{ref_id}", authors:{author_string})\n'

    return s


def create_reference_mutation(ref_id, ref):
    s = f'''{ref_id}: createLiteratureReference(id: "{ref_id}", title: "{ref['title']}", volume: "{ref['volume']}", first_page: "{ref['first_page']}", last_page: "{ref['last_page']}", publication_Year: "{ref['year']}", DOI: "{ref['DOI']}", PMID: "{ref['PubMed']}")
'''
    return s


def create_author_mutation(id,surname,first,middle):
    s = f'''{id}: createAuthor(id: "{id}",surname: "{surname}",first_initial: "{first}" ,middle_initial: "{middle}")
'''
    return s


def create_journal_mutation(journal, journal_id):
    s = f'''{journal_id}: createJournal(id: "{journal_id}",name: "{journal}")
'''
    return s


def write_mutation(io_list, pheno_list):
    pheno_dict = {}
    biomarker_dict = {}
    author_dict = {}
    journal_dict = {}
    reference_dict = {}
    location_dict = {}

    muatation_number = 0
    s = get_mutation_start(muatation_number)
    s += write_phenotype_mutation(pheno_list,pheno_dict)
    s += '}\n'
    print(s)
    muatation_number += 1

    for index, io in enumerate(io_list, start=1):
        s = get_mutation_start(muatation_number)
        gene_id,m = write_gene_mutation(index, io)
        s += m
        isoform,m = write_isoform_mutation(index, io,gene_id)
        s += m
        protein,m = write_protein_mutation(index, io,gene_id,isoform)
        s += m
        if io['Marker'] in biomarker_dict:
            biomarker = biomarker_dict[io['Marker']]
        else:
            biomarker,m = write_biomarker_mutation(index,io)
            s += m
            biomarker_dict[io['Marker']] = biomarker
        m = write_biomarker_protein_mutation(biomarker,protein)
        s += m
        phenotype = pheno_dict[io['ImmunePhenotype']]
        m = write_biomarker_phenotype_mutation(biomarker, phenotype)
        s += m
        m = write_location_mutations(io['sp_info'],author_dict,journal_dict,reference_dict,location_dict,biomarker)
        s += m
        s += '}\n'
        print(s)
        muatation_number += 1


def get_mutation_start(muatation_number):
    mutation_id = f'm{muatation_number}'
    s = 'mutation ' + str(mutation_id) + ' {\n'
    return s


def write_biomarker_mutation(index, io):
    s = ''
    id = 'b' + str(index)
    biomarker = io['Marker']
    s = f'{id}: createBioMarker(id: "{id}", name: "{biomarker}", synonyms: [])\n'
    return id,s



def write_biomarker_protein_mutation(biomarker, protein):
    s = ''
    # AddBioMarkerPolypeptide_chains
    tag = biomarker + '_' + protein
    s = f'{tag}: addBioMarkerPolypeptide_chains(id:"{biomarker}", polypeptide_chains:["{protein}"])\n'
    return s


# AddImmunePhenotypeBiomarkers(
# from: _ImmunePhenotypeInput!
# to: _BioMarkerInput!
# ):
def write_biomarker_phenotype_mutation(biomarker, phenotype):
    # AddBioMarkerPolypeptide_chains
    tag = biomarker + '_' + phenotype
    s = f'{tag}: addImmunePhenotypeBiomarkers(id:"{phenotype}", biomarkers:["{biomarker}"])\n'
    return s


def write_protein_mutation(index, io,gene_id,isoform):
    sp_info = io['sp_info']
    id = 'p' + str(index)
    protein = sp_info['name']
    uniprot_id = sp_info['id']
    acc_num = sp_info['acc_num']
    s = f'{id}: createProtein(id: "{id}", name: "{protein}", uniprot_id: "{uniprot_id}", accessionNumber: "{acc_num}" )\n'

    # AddProteinGene
    tag = id + '_' + gene_id
    s += f'{tag}: addProteinGene(gene:"{gene_id}", id:"{id}")\n'

    # AddProteinTranscript
    tag = id + '_' + isoform
    s += f'{tag}: addProteinTranscript(transcript:"{isoform}", id:"{id}")\n'
    return id,s


def write_isoform_mutation(index, io,gene_id):
    id = 'i' + str(index)
    isoform = io['sp_info']['transcript']
    s = f'{id}: createIsoform(id: "{id}", name: "{isoform}", )\n'

    tag = id + '_' + gene_id

    s += f'{tag}: addIsoformGene(gene:"{gene_id}",id:"{id}")\n'
    return id, s


# CreateGene(
# id: ID
# name: String!
# chromosome: String!
# strand: Strand!
# start: Int
# end: Int

def write_gene_mutation(index, io):
    id = 'g' + str(index)
    gene_info = io['gene_info']
    gene = io['GeneSymbol']
    chrom = gene_info['chrom']
    strand = gene_info['strand']
    start = gene_info['start']
    end = gene_info['end']
    s = f'{id}: createGene(id: "{id}", name: "{gene}", chromosome: "{chrom}", strand:{strand}, start: {start}, end: {end})\n'
    ref = gene_info['reference']
    if ref['type']=='InternetReference':
        ref_id = 'ref_' + id
        s += f'{ref_id}: createInternetReference(id:"{ref_id}", web_address:"{ref["url"]}", accessed_date:"02-22_2020")\n'
        ref_id2 = 'gref_' + id
        s += f'{ref_id2}: addGeneReference(id:"{id}", reference:"{ref_id}" )\n'
    return id,s


def write_phenotype_mutation(pheno_list, pheno_dict):
    s = ''
    for index, pheno in enumerate(pheno_list, start=1):
        id = 'ph' + str(index)
        pheno_dict[pheno] = id
        s += f'{id}: createImmunePhenotype(id: "{id}", name: "{pheno}")\n'
    return s



def main():
    io_list, pheno_list = get_io_and_pheno_lists()

    write_mutation(io_list,pheno_list)
    # out_file = open('data/json_out', 'w+')
    # json.dump(io_list, out_file)


    # with open('data/json_out') as fh:
    #     a = json.load(fh)
    #     print(a)

    # write_neo_graph(io_list, pheno_list)


if __name__ == "__main__":
    main()