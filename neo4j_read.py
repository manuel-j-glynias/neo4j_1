
import json

from neo4j import GraphDatabase
import csv
import uniprot

def get_io_list():
    # file_path = 'data/Small_MarkerReporting.csv'
    file_path = 'data/tblPROD_OA_OA_MarkerReporting.csv'
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


# MATCH (a:Biomarker)-[:CodedBy]->(b)
#  WHERE b.name = 'LAG3'
#  RETURN a,b
def print_coded_by(tx, name):
    for record in tx.run("MATCH (a:Biomarker)-[:CodedBy]->(g) WHERE g.name = {name} RETURN a.name", name=name):
        print(record["a.name"])



def create_phenotype(tx, pheno, id):
    tx.run("CREATE (p:ImmunePhenotype{id:{id}, name:{phenotype}})",id=id,phenotype=pheno)


def create_gene(tx, gene, id,chromosome,strand,start,end):
    tx.run("CREATE (g:Gene{id:{id}, name:{gene}, chromosome: {chromosome}, strand: {strand}, start: {start}, end: {end}})",
                  id=id,gene=gene,chromosome=chromosome,strand=strand,start=start,end=end)


def create_isoform(tx, name, id,gene):
    tx.run("CREATE (i:Isoform{id:{id}, name:{name}})",id=id,name=name)
    tx.run("MATCH (i:Isoform{name:{name}}) MATCH (g:Gene{name:{gene}}) MERGE (i)<-[:TRASCRIBED_FROM]-(g)",name=name,gene=gene)


def create_protein(tx, name, id,uniprot_id,accnum,isoform,gene):
    tx.run("CREATE (i:Protein{id:{id}, name:{name}, uniprot_id:{uniprot_id}, accessionNumber:{accnum} })",
           id=id,name=name,uniprot_id=uniprot_id,accnum=accnum)
    tx.run("MATCH (p:Protein{name:{name}}) MATCH (g:Gene{name:{gene}}) MERGE (p)<-[:CODED_BY]-(g)",name=name,gene=gene)
    tx.run("MATCH (p:Protein{name:{name}}) MATCH (i:Isoform{name:{isoform}}) MERGE (p)<-[:TRANSLATED_FROM]-(i)",name=name,isoform=isoform)


def create_biomarker(tx,biomarker, id):
    tx.run("CREATE (b:BioMarker{id:{id}, name:{name}})",id=id,name=biomarker)

def biomarker_exists(tx, name):
    b = False
    for record in tx.run("MATCH (a:BioMarker) WHERE a.name = {name} RETURN a.name", name=name):
        b = True
    return b

def relate_biomarker(tx,biomarker,phenotype,protein):
    tx.run("MATCH (b:BioMarker{name:{name}}) MATCH (p:Protein{name:{protein}}) MERGE (b)<-[:MONOMER]-(p)",name=biomarker,protein=protein)
    tx.run("MATCH (ph:ImmunePhenotype{name:{phenotype}}) MATCH (b:BioMarker{name:{biomarker}}) MERGE (b)-[:PHENOTYPE]->(ph)",
           biomarker=biomarker,phenotype=phenotype)


def get_io_and_pheno_lists():
    io_list, pheno_list = get_io_list()
    for io in io_list:
        gene_info, sp_info = uniprot.add_gene_and_sp_info(io['GeneSymbol'])
        io['gene_info'] = gene_info
        io['sp_info'] = sp_info
        print(io['GeneSymbol'], io['sp_info']['transcript'])
    return io_list, pheno_list


def write_neo_graph(io_list, pheno_list):
    uri = "bolt://localhost:7687"
    driver = GraphDatabase.driver(uri, auth=("neo4j", "omni"))
    with driver.session() as session:
        tx = session.begin_transaction()
        for index, pheno in enumerate(pheno_list, start=1):
            id = 'ph' + str(index)
            create_phenotype(tx, pheno, id)
        for index, io in enumerate(io_list, start=1):
            id = 'g' + str(index)
            gene_info = io['gene_info']
            gene = io['GeneSymbol']
            create_gene(tx, gene, id, gene_info['chrom'], gene_info['strand'], gene_info['start'], gene_info['end'])
            id = 'i' + str(index)
            isoform = io['sp_info']['transcript']
            create_isoform(tx, isoform, id, gene)
            sp_info = io['sp_info']
            id = 'p' + str(index)
            create_protein(tx, sp_info['name'], id, sp_info['id'], sp_info['acc_num'], isoform, gene)
            if not biomarker_exists(tx, io['Marker']):
                id = 'b' + str(index)
                create_biomarker(tx, io['Marker'], id)
            relate_biomarker(tx, io['Marker'], io['ImmunePhenotype'], sp_info['name'])

        tx.commit()

def write_mutation(io_list, pheno_list):
    pheno_dict = {}
    s = """mutation {
        """
    s += write_phenotype_mutation(pheno_list,pheno_dict)
    biomarker_dict = {}
    for index, io in enumerate(io_list, start=1):
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
    s += '''
    }'''
    print(s)


def write_biomarker_mutation(index, io):
    s = ''
    id = 'b' + str(index)
    biomarker = io['Marker']
    s = f'''
    {id}: CreateBioMarker(
        id: "{id}"
        name: "{biomarker}"
        synonyms: "" 
        )'''
    s += """
        {
           id
           name
         }
    """
    return id,s



def write_biomarker_protein_mutation(biomarker, protein):
    s = ''
    # AddBioMarkerPolypeptide_chains
    tag = biomarker + '_' + protein
    s += tag
    s += ''': AddBioMarkerPolypeptide_chains(
         from: { 
             id:'''
    s += f'"{biomarker}"'
    s += '''
             }
         to: { 
             id:'''
    s += f'"{protein}"'

    s += """
             }
         )
             {
              from {
                 id
                 name
             }
             to {
                 id
                 name
             }
           }
           """
    return s
# AddImmunePhenotypeBiomarkers(
# from: _ImmunePhenotypeInput!
# to: _BioMarkerInput!
# ):
def write_biomarker_phenotype_mutation(biomarker, phenotype):
    s = ''
    # AddBioMarkerPolypeptide_chains
    tag = biomarker + '_' + phenotype
    s += tag
    s += ''': AddImmunePhenotypeBiomarkers(
         from: { 
             id:'''
    s += f'"{phenotype}"'
    s += '''
             }
         to: { 
             id:'''
    s += f'"{biomarker}"'

    s += """
             }
         )
             {
              from {
                 id
                 name
             }
             to {
                 id
                 name
             }
           }
           """
    return s


def write_protein_mutation(index, io,gene_id,isoform):
    sp_info = io['sp_info']
    id = 'p' + str(index)
    protein = sp_info['name']
    uniprot_id = sp_info['id']
    acc_num = sp_info['acc_num']
    s = f'''
    {id}: CreateProtein(
        id: "{id}"
        name: "{protein}"
        uniprot_id: "{uniprot_id}"
        accessionNumber: "{acc_num}" 
        )'''
    s += """
        {
           id
           name
         }
         """
    # AddProteinGene
    tag = id + '_' + gene_id
    s += tag
    s += ''': AddProteinGene(
        from: { 
            id:'''
    s += f'"{id}"'
    s += '''
            }
        to: { 
            id:'''
    s += f'"{gene_id}"'

    s += """
            }
        )
            {
             from {
                id
                name
            }
            to {
                id
                name
            }
          }
          """
    # AddProteinTranscript
    tag = id + '_' + isoform
    s += tag
    s += ''': AddProteinTranscript(
        from: { 
            id:'''
    s += f'"{id}"'
    s += '''
            }
        to: { 
            id:'''
    s += f'"{isoform}"'

    s += """
            }
        )
            {
             from {
                id
                name
            }
            to {
                id
                name
            }
          }
          """

    return id,s


def write_isoform_mutation(index, io,gene_id):
    id = 'i' + str(index)
    isoform = io['sp_info']['transcript']
    s = f'''
    {id}: CreateIsoform(
        id: "{id}"
        name: "{isoform}"
        )'''
    s += """
        {
           id
           name
         }
    """
    tag = id + '_' + gene_id
    s += tag
    s += ''': AddIsoformGene(
        from: { 
            id:'''
    s += f'"{id}"'
    s += '''
            }
        to: { 
            id:'''
    s += f'"{gene_id}"'

    s += """
            }
        )
            {
             from {
                id
                name
            }
            to {
                id
                name
            }
          }
          """
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
    s = f'''
    {id}: CreateGene(
        id: "{id}"
        name: "{gene}"
        chromosome: "{chrom}"
        strand:{strand}
        start: {start}
        end: {end}
        )'''
    s += """
        {
           id
           name
         }
         """
    return id,s


def write_phenotype_mutation(pheno_list, pheno_dict):
    s = ''
    for index, pheno in enumerate(pheno_list, start=1):
        id = 'ph' + str(index)
        pheno_dict[pheno] = id
        s += f'''
    {id}: CreateImmunePhenotype(
        id: "{id}"
        name: "{pheno}"
        )'''
        s += """
        {
            id
            name
        }"""
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