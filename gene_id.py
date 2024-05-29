#!/usr/bin/env python3
# Author: Jiwoong Kim (jiwoongbio@gmail.com)

import re
import os
import sys
import time
import gzip
import argparse
import requests
import warnings
import xml.etree.ElementTree as ET
import urllib.parse
from openai import OpenAI

def download_file(url, file_path):
	with warnings.catch_warnings():
		warnings.simplefilter("ignore", requests.urllib3.exceptions.InsecureRequestWarning)
		response = requests.get(url, stream = True, verify = False)
		response.raise_for_status()
		with open(file_path, 'wb') as file:
			for chunk in response.iter_content(chunk_size = 8192):
				if chunk:
					file.write(chunk)

def download_json(url):
	headers = {
		'Accept': 'application/json',
	}
	with warnings.catch_warnings():
		warnings.simplefilter("ignore", requests.urllib3.exceptions.InsecureRequestWarning)
		response = requests.get(url, headers = headers, verify = False)
		response.raise_for_status()
		return response.json()

def download_content(url):
	with warnings.catch_warnings():
		warnings.simplefilter("ignore", requests.urllib3.exceptions.InsecureRequestWarning)
		response = requests.get(url, verify = False)
		response.raise_for_status()
		return response.content

def get_child_node_list(node, tag_list):
	child_node_list = []
	tag = tag_list.pop(0)
	for child_node in node:
		if child_node.tag == tag:
			if tag_list:
				child_node_list.extend(get_child_node_list(child_node, tag_list))
			else:
				child_node_list.append(child_node)
	return child_node_list

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("gene_query", nargs = "*", help = "gene queries")
	parser.add_argument("-r", "--redownload", dest = "redownload", help = "redownload data", action = "store_true")
	parser.add_argument("-t", "--taxonomy_id", dest = "taxonomy_id", help = "NCBI taxonomy ID", type = int)
	args = parser.parse_args()

	directory = os.path.dirname(os.path.abspath(__file__))
	gene_info_path = os.path.join(directory, "gene_info.gz")
	gene_id_symbol_path = os.path.join(directory, "gene_id.symbol.txt")
	gene_id_synonym_path = os.path.join(directory, "gene_id.synonym.txt")

	if args.redownload or not os.path.exists(gene_info_path):
		download_file("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz", gene_info_path)

	if args.redownload or not os.path.exists(gene_id_symbol_path):
		if args.taxonomy_id is not None:
			with gzip.open(os.path.join(directory, "gene_info.gz"), 'rt') as input_file, open(gene_id_symbol_path, 'w') as output_file:
				columns = re.sub(r'^#', '', input_file.readline()).strip().split('\t')
				for line in input_file:
					line_dict = dict(zip(columns, line.strip().split('\t')))
					if int(line_dict['tax_id']) == args.taxonomy_id:
						output_file.write('\t'.join([line_dict['GeneID'], line_dict['Symbol']]) + '\n')

	symbol_gene_id_list_dict = {}
	with open(gene_id_symbol_path, 'r') as input_file:
		for line in input_file:
			(gene_id, symbol) = line.strip().split('\t')
			gene_id = int(gene_id)
			symbol_gene_id_list_dict.setdefault(symbol, []).append(gene_id)

	if args.redownload or not os.path.exists(gene_id_synonym_path):
		if args.taxonomy_id is not None:
			with gzip.open(os.path.join(directory, "gene_info.gz"), 'rt') as input_file, open(gene_id_synonym_path, 'w') as output_file:
				columns = re.sub(r'^#', '', input_file.readline()).strip().split('\t')
				for line in input_file:
					line_dict = dict(zip(columns, line.strip().split('\t')))
					if int(line_dict['tax_id']) == args.taxonomy_id:
						for symbol in line_dict['Synonyms'].split('|') + line_dict['LocusTag'].split('|'):
							if symbol in symbol_gene_id_list_dict:
								continue
							if symbol == '-':
								continue
							output_file.write('\t'.join([line_dict['GeneID'], symbol]) + '\n')

	with open(gene_id_synonym_path, 'r') as input_file:
		for line in input_file:
			(gene_id, symbol) = line.strip().split('\t')
			gene_id = int(gene_id)
			symbol_gene_id_list_dict.setdefault(symbol, []).append(gene_id)

	gene_query_list = []
	if args.gene_query is not None:
		for gene_query in args.gene_query:
			if gene_query == '-':
				for line in sys.stdin:
					gene_query_list.append(line.strip())
			else:
				gene_query_list.append(gene_query)

	for gene_query in gene_query_list:
		symbol_list = re.split(r', *', gene_query)
		if any(' ' in symbol for symbol in symbol_list):
			symbol_list = extract_symbol_list(gene_query)
		for symbol in symbol_list:
			gene_id_list = []
			if symbol == '':
				pass
			elif symbol in symbol_gene_id_list_dict:
				gene_id_list.extend(symbol_gene_id_list_dict[symbol])
			elif symbol.startswith("ncbigene:"):
				gene_id = symbol.split(":", 1)[1]
				gene_id_list.append(gene_id)
			elif symbol.startswith("hgnc:"):
				hgnc_id = symbol.split(":", 1)[1]
				url = f"https://rest.genenames.org/fetch/hgnc_id/{hgnc_id}"
				data = download_json(url)
				gene_id_list.extend([int(x['entrez_id']) for x in data['response']['docs']])
			elif symbol.startswith("ensembl:"):
				ensembl_id = symbol.split(":", 1)[1]
				url = f"https://rest.ensembl.org/xrefs/id/{ensembl_id}?external_db=EntrezGene;content-type=application/json"
				data = download_json(url)
				gene_id_list.extend([int(x['primary_id']) for x in data if 'primary_id' in x])
			else:
				if args.taxonomy_id is not None:
					term = f"\"{symbol}\"[Gene Name] AND txid{args.taxonomy_id}[Organism]"
					gene_id_list.extend(esearch("gene", term))
			print('\t'.join([gene_query, symbol, ','.join([str(gene_id) for gene_id in gene_id_list])]))

def esearch(db, term):
	api_key = os.environ.get("NCBI_API_KEY")
	term = urllib.parse.quote(term)
	id_list = []
	if api_key is not None:
		url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={term}&api_key={api_key}"
		time.sleep(0.105)
	else:
		url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={db}&term={term}"
		time.sleep(0.350)
	root = ET.fromstring(download_content(url))
	id_list.extend([int(x.text) for x in get_child_node_list(root, ["IdList", "Id"])])
	return id_list

def extract_symbol_list(gene_query):
	client = OpenAI(
		api_key = os.environ.get("OPENAI_API_KEY"),
	)
	completion = client.chat.completions.create(
		model = "gpt-4o",
		messages = [
			{"role": "user", "content": f"Extract all possible gene names from \"{gene_query}\" and absolutely output them as a comma-separated list only. Do not guess. If there are no gene names, absolutely output 0 only."},
		],
		temperature = 0.2,
	)
	gene_query = completion.choices[0].message.content
	if gene_query == "0":
		gene_query = ''
	symbol_list = re.split(r', *', gene_query)
	return symbol_list

if __name__ == "__main__":
	main()
