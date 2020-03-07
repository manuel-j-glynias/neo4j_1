import requests, sys
import datetime

def fetch_uniprot_by_gene_name(gene_name):
  requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&taxid=9606&exact_gene=" + gene_name
  r = requests.get(requestURL, headers={ "Accept" : "application/json"})
  if not r.ok:
    r.raise_for_status()
    sys.exit()

  responseBody = r.json()
  # pprint.pprint(responseBody)
  return responseBody[0]


def fetch_uniprot_by_acc_num(id):
  requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession=" + id
  r = requests.get(requestURL, headers={ "Accept" : "application/json"})
  if not r.ok:
    r.raise_for_status()
    sys.exit()
  responseBody = r.json()
  # pprint.pprint(responseBody)
  return responseBody[0]


def get_from__dbReferences(ref_list,type):
  val = None
  for ref in ref_list:
      if ref['type'] == type:
        val = ref['id']
        break
  return val


def get_from__comments(comments_list,type,identifier):
  val = None
  for comment in comments_list:
      if comment['type'] == type:
        val = ''
        for c in comment[identifier]:
          val += c['value']
        break
  return val


def get_location_references_from_uniport_id(uniprot, ref_list):
  sp = fetch_uniprot_by_acc_num(uniprot)
  references = get_reference_list_for_scope(sp['references'], 'SUBCELLULAR LOCATION')
  ref_list.extend(references)


def get_reference_list_for_scope(ref_list, scope):
  reference_list = []
  for ref in ref_list:
    if scope in ref['scope']:
      cit = ref['citation']
      if cit['type'] == 'journal article':
        reference = extract_reference(cit)
        reference['evidence_type'] = 'sequence similarity evidence used in manual assertion'
        reference_list.append(reference)
  return reference_list


def get_location_from__comments(comments_list,reference_dict,acc_num):
  val = None
  ref_list = []
  for comment in comments_list:
      if comment['type'] == 'SUBCELLULAR_LOCATION':
        val = ''
        for c in comment['locations']:
          val += c['location']['value']
          if 'evidences' in c['location']:
            for evidence in c['location']['evidences']:
              if evidence['code']=='ECO:0000269' and evidence['source']['name']== 'PubMed':
                pubmed = evidence['source']['id']
                if pubmed in reference_dict:
                  ref = reference_dict[pubmed]
                  ref['evidence_type'] = 'experimental evidence used in manual assertion'
                  ref_list.append(ref)
                # else:
                #   print("miss")
              elif evidence['code']=='ECO:0000250':
                if 'source' in evidence and evidence['source']['name']== 'UniProtKB':
                  uniprot = evidence['source']['id']
                  get_location_references_from_uniport_id(uniprot,ref_list)
              # elif evidence['code']=='ECO:0000305':
              #   print('305')
  if val is not None and len(ref_list)==0:
    ref_list.append(create_uniprot_reference(acc_num))
  return val,ref_list


def get_transcript_from__dbReferences(ref_list):
  val = None
  for ref in ref_list:
      if ref['type'] == 'RefSeq':
        if 'properties' in ref:
          val = ref['properties']['nucleotide sequence ID']
        break
  return val


def get_pathways_from__dbReferences(ref_list):
  pathways = []
  for ref in ref_list:
      if ref['type'] == 'Reactome':
        if 'properties' in ref:
          name = ref['properties']['pathway name']
          pathways.append(name)
  return pathways


def get_reference_dict_for_scope(ref_list, scope):
  reference_dict = {}
  for ref in ref_list:
    if scope in ref['scope']:
      cit = ref['citation']
      if cit['type'] == 'journal article':
        reference = extract_reference(cit)
        reference_dict[reference['PubMed']] = reference
  return reference_dict


def extract_reference(cit):
  reference = {'DOI':'_'}
  reference['type'] = 'LiteratureReference'
  reference['title'] = cit['title']
  reference['authors'] = cit['authors']
  reference['journal'] = cit['publication']['journalName']
  reference['year'] = cit['publicationDate']
  reference['volume'] = cit['location']['volume']
  reference['first_page'] = cit['location']['firstPage']
  reference['last_page'] = cit['location']['lastPage']
  for dbRef in cit['dbReferences']:
    if dbRef['type'] == 'PubMed':
      reference['PubMed'] = dbRef['id']
    elif dbRef['type'] == 'DOI':
      reference['DOI'] = dbRef['id']
  # if reference['PubMed']=='28484017':
  #   print('here')
  return reference


def create_uniprot_reference(acc_num):
  reference = {}
  reference['type'] = 'InternetReference'
  reference['url'] = 'https://www.uniprot.org/uniprot/' + acc_num
  reference['accessed_date'] = datetime.datetime.now()
  return reference


def create_mygene_reference(gene_id):
  reference = {}
  reference['type'] = 'InternetReference'
  reference['url'] = "http://mygene.info/v3/gene/" + gene_id
  reference['accessed_date'] = datetime.datetime.now()
  return reference


def get_sp_info(id,gene_name):
  sp_info = {'function':'', 'location':'', 'transcript':''}
  if len(id)>0:
    sp = fetch_uniprot_by_acc_num(id)
  else:
    sp = fetch_uniprot_by_gene_name(gene_name)
  sp_info['acc_num'] = sp['accession']
  sp_info['id'] = sp['id']
  reference_dict = get_reference_dict_for_scope(sp['references'], 'SUBCELLULAR LOCATION')

  if 'recommendedName' in sp['protein']:
    sp_info['name'] = sp['protein']['recommendedName']['fullName']['value']
  elif 'submittedName' in sp['protein']:
    sp_info['name'] = sp['protein']['submittedName'][0]['fullName']['value']
  if 'comments' in sp:
    sp_info['function'] = get_from__comments(sp['comments'],'FUNCTION','text')
    sp_info['location'],sp_info['location_references'] = get_location_from__comments(sp['comments'],reference_dict,sp_info['acc_num'])
  sp_info['transcript'] = get_transcript_from__dbReferences(sp['dbReferences'])
  sp_info['pathways'] = get_pathways_from__dbReferences(sp['dbReferences'])

  # print('location=',sp_info['location'])
  # if len(sp_info['location_references'])==0:
  #   print("none")
  # print(sp_info['location_references'])

  # sp_info['geneID'] = get_from__dbReferences(sp['dbReferences'], 'GeneID')
  return sp_info


def fetch_gene_id_by_gene_name(gene_name):
  requestURL = "http://mygene.info/v3/query?species=human&q=symbol:" + gene_name
  r = requests.get(requestURL, headers={ "Accept" : "application/json"})
  if not r.ok:
    r.raise_for_status()
    sys.exit()

  responseBody = r.json()
  # pprint.pprint(responseBody)
  return responseBody['hits'][0]['entrezgene']

def fetch_gene_info_by_gene_id(gene_id):
  requestURL = "http://mygene.info/v3/gene/" + gene_id
  r = requests.get(requestURL, headers={ "Accept" : "application/json"})
  if not r.ok:
    r.raise_for_status()
    sys.exit()

  responseBody = r.json()
  # pprint.pprint(responseBody)
  return responseBody


def get_gene_info(gene_name):
  gene_id = fetch_gene_id_by_gene_name(gene_name)
  gene_info = fetch_gene_info_by_gene_id(gene_id)
  gene_dict = {'summary':'', 'sp':''}
  genomic_pos = gene_info['genomic_pos_hg19']
  if isinstance(genomic_pos, list):
    genomic_pos = genomic_pos[0]
  gene_dict['chrom'] = genomic_pos['chr']
  gene_dict['strand'] = 'FORWARD'
  if genomic_pos['strand']== -1:
    gene_dict['strand'] = 'REVERSE'
  gene_dict['start'] = genomic_pos['start']
  gene_dict['end'] = genomic_pos['end']
  gene_dict['reference'] = create_mygene_reference(gene_id)
  if 'summary' in gene_info:
    gene_dict['summary'] = gene_info['summary']
  if 'uniprot' in gene_info:
    if 'Swiss-Prot' in gene_info['uniprot']:
      gene_dict['sp'] = gene_info['uniprot']['Swiss-Prot']
  return gene_dict


def add_gene_and_sp_info(target, gene_dict, protein_dict):
  gene_symbol_dict = {'CD8':['CD8A','CD8B']}
  target['components'] = []
  gene_array = [target['GeneSymbol']]
  if target['GeneSymbol'] in gene_symbol_dict:
    gene_array = gene_symbol_dict[target['GeneSymbol']]
  for gene_name in gene_array:
    component = {}
    if gene_name in gene_dict:
      gene_info = gene_dict[gene_name]
    else:
      gene_info = get_gene_info(gene_name)
      gene_dict[gene_name] = gene_info
    if gene_info['sp'] in protein_dict:
      sp_info = protein_dict[gene_info['sp']]
    else:
      sp_info = get_sp_info(gene_info['sp'], gene_name)
      protein_dict[gene_info['sp']] = sp_info

    component['GeneSymbol'] = gene_name
    component['gene_info'] = gene_info
    component['sp_info'] = sp_info
    target['components'].append(component)


def main():
  gene_info, sp_info = add_gene_and_sp_info('CD4')
  print(sp_info)
  print(gene_info)


if __name__ == "__main__":
    main()