import argparse
import os, sys
import json
import wget
import shutil 
import gzip
import subprocess
import pandas as pd
import multiprocessing as mp
from os import path, listdir
from collections import Counter
os.chmod('./scripts/liftOver', 0o755)
### Options from command line
parser = argparse.ArgumentParser(description = "ConnectOR (v.2.0)"\
                                               "By: Carlos Pulido (carlos.pulidoquetglas@unil.ch)")

def range_limited_int_type(arg):
    """ Type function for argparse - an int within some predefined bounds """
    try:
        i = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer")
    if i < 1 or i > 99:
        raise argparse.ArgumentTypeError("Argument must be < " + str(1) + " and > " + str(99) + ". Default: 30")
    return i

parser.add_argument('-mM', '--minMatch',
                    required = False,
                    default = '30',
                    action = 'store',
                    dest = 'minMatch',
                    type = range_limited_int_type,
                    help = '0.N Minimum ratio of bases that must remap in liftOver step.'\
                           'Default: 30 (0.30 minimum ratio)')

parser.add_argument('-g', '--gene',
                    dest = 'gene_level',
                    action = 'store_true',
                    default = False,
                    help = 'Generate results at gene level along with exon level results (default: False)')

parser.add_argument('-r', '--redo',
                    dest = 'redo',
                    action = 'store_true',
                    default = False,
                    help = 'Force pipeline to redo all steps in the analysis (default: False)')


### Storing user options
options = parser.parse_args()
gene_level = options.gene_level
minMatch = options.minMatch


### Dictionaris
with open('dictionaries.json') as f:
    dictionaries = json.load(f)


### FUNCTIONS
def assign_simplified_geneType(Type):
    type_dict = {"none": ["none"],
                 "pc": ["protein_coding", "pc"],

                "lncRNA": ["lncRNA", "TEC", "processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",
"non_coding", "sense_intronic", "sense_overlapping", "macro_lncRNA","antisense_RNA", 
"bidirectional_promoter_lncrna", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA"],

                 "ncRNA": ["ncRNA", "snRNA", "snoRNA", "rRNA", "rRNA_pseudogene", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "miRNA", "ribozyme", "sRNA", 
"scaRNA", "scRNA", "vault_RNA", "known_ncrna","vaultRNA"],

                 "other": ["other", "processed_pseudogene", "transcribed_processed_pseudogene", "translated_processed_pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "transcribed_unitary_pseudogene", "polymorphic pseudogene", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "IG_C_pseudogene", "IG_D_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "IG_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "NOVEL"]}
    biotype = [t for t in type_dict if Type in type_dict[t]]
    #if biotype == []: biotype = ["ncRNA"]
    if biotype == []: biotype = ["other"]
    return(biotype[0])


def read_config(file_name):
    config_df =  pd.read_csv(file_name, sep='\t', na_filter= False)
    return(config_df)


def check_folder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


def check_file(file_name_list):
    for file_name in file_name_list:
        if not os.path.isfile(file_name):
            print("'%s' doesn not exist. Please check config_file. Exiting..."%(file_name))
            sys.exit(1)


def change_default(feature, value, config_df, i):
    default_values = config_df.loc[i]["default"].split("|")
    if feature == "gtf":
        n = 0
    else:
        n = 1
    default_values[n] = value
    return("|".join(default_values))


def download_default_files(file_url, folder_name):
    check_folder(folder_name)
    file_name = file_url.split("/")[-1]
    if not os.path.isfile("/".join([folder_name, file_name])):
        wget.download(file_url, folder_name)
        print('\t%s downloaded succesfully..'%(file_name))
    else:
        print('\t%s already exists.. Skipping'%(file_name))


def arguments_dic(line):
    args = line[-1].split(";")
    dic = {}
    for e in args:
        if e == ";": continue
        if e == "": continue
        if e[0] == " ": e = e[1:]
        key = e.split(" ")[0].replace("\"","")
        dic[key] = e.split(" ")[1].replace("\"","")
    return dic


def biotype(keys, arguments):
    biotype=''
    while not biotype:
        for k in keys:
            try:
                biotype = arguments[k]
                return(biotype)
            except KeyError:
                pass
        if not biotype: biotype = "NOVEL"
        return(biotype)


#def generate_maps(gtf, sp):
def generate_maps(arg_list):
    gtf = arg_list[0]
    sp = arg_list[1]
    if not options.redo and os.path.isfile("maps/"+sp+".geneID_map.txt"):
        print("\tMaps for %s already exist.. Skipping"%(sp))
        return
    else:
        print("\tMaps for %s... generating"%(sp))

    if gtf.endswith(".gz"):
        f = gzip.open(gtf, 'rb')
        compressed = True
    else:
        f = open(gtf, 'r')
        compressed = False

    genes = {}
    for line in f:
        if compressed: line = str(line, 'utf-8')
        if line.startswith("#"): continue
        line = line.strip().split("\t")
        arguments = arguments_dic(line)

        gene_biotype = biotype(["gene_type", "gene_biotype"], arguments)
        gene_biotype = assign_simplified_geneType(gene_biotype)

        if "gene_name" in arguments:
            gene_name = arguments["gene_name"]
            if gene_name.startswith("lncRNA:" or "lincRNA:"):
                gene_biotype = 'lncRNA'
        else:
            gene_name = arguments["gene_id"]

        geneID = arguments["gene_id"]

        if line[2] == "exon":
            if not geneID in genes:
                genes[geneID] = {"gene_name": gene_name,
                                "gene_type": gene_biotype,
                                "coordinates": [line[0],line[3],line[4],line[6]],
                                "transcripts": []}

            transcriptID = arguments["transcript_id"]
            if not transcriptID in genes[geneID]["transcripts"]:
                genes[geneID]["transcripts"].append(transcriptID)

            exon_start = line[3]
            exon_end = line[4]
            if int(exon_start) < int(genes[geneID]["coordinates"][1]): genes[geneID]["coordinates"][1] = exon_start
            if int(exon_end) > int(genes[geneID]["coordinates"][2]): genes[geneID]["coordinates"][2] = exon_end

            # Update gene biotype; protein_coding prevails other biotypes
            if genes[geneID]["gene_type"] == gene_biotype or genes[geneID]["gene_type"] == "protein_coding":
                continue
            else:
                genes[geneID]["gene_type"] = gene_biotype

    fo = open("maps/"+sp+".geneID_map.txt", 'w')
    for g in genes:
        line = "\t".join([g, genes[g]["gene_name"], genes[g]["gene_type"], ":".join(genes[g]["coordinates"]), "|".join(genes[g]["transcripts"])])
        fo.write(line+"\n")
    fo.close()
    print("\tMaps for %s... generated "%(sp))


def generate_beds(arg_list):
    gtf = arg_list[0]
    sp = arg_list[1]
    gene_level = arg_list[2]
    
    if gtf.endswith(".gz"):
        f = gzip.open(gtf, 'rb')
        compressed = True
    else:
        f = open(gtf, 'r')
        compressed = False

    # Check if exon BEDs are required
    if not options.redo and os.path.isfile("./BEDs/{}.exons.bed".format(sp)):
        print("\tBED for {} exons already exists.. Skipping".format(sp))
        redo_exons = False
    else:
        print("\tExons BEDs for {}... generating".format(sp))
        output_exons = open("./BEDs/{}.exons.bed".format(sp), "w")
        redo_exons = True

    # Check if genes BEDs are required
    if not options.redo and os.path.isfile("./BEDs/{}.genes.bed".format(sp)):
        print("\tBED for {} genes already exists.. Skipping".format(sp))
        redo_genes = False
    else:
        print("\tGenes BEDs for {}... generating".format(sp))
        output_genes = open("./BEDs/{}.genes.bed".format(sp), "w")
        redo_genes = True

    # Skip function if outputs already exist
    if not redo_exons and not redo_genes:
        return

    # Read GTF and create output BEDs
    genes = {}
    for line in f:
        if compressed:
            line = str(line, 'utf-8')
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        if line[2] != "exon":
            continue
        arguments = arguments_dic(line)

        if "gene_name" in arguments:
            gene_name = arguments["gene_name"]
        else:
            gene_name = arguments["gene_id"]

        gene_id = arguments["gene_id"]
        gene_biotype = arguments.get("gene_type") or arguments.get("gene_biotype")

        # Handle gene biotype assignment
        if gene_name.startswith(("lncRNA:", "lincRNA:")):
            gene_biotype = 'lncRNA'  # Assign default value if not found

        gene_biotype = assign_simplified_geneType(gene_biotype)

        chrom = line[0] if line[0].startswith("chr") else "chr" + line[0]
        start = str(int(line[3]) - 1)
        end = line[4]
        strand = line[6]

        if int(start) < 0:
            start = "0"
        coordinates = ":".join([chrom, start, end, strand])
        exon_id = "|".join([gene_id, gene_name, gene_biotype, coordinates])

        if redo_exons:
            exon_bed_line = "\t".join([chrom, start, end, exon_id, '0', strand]) + "\n"
            output_exons.write(exon_bed_line)

        if redo_genes:
            if not gene_id in genes:
                genes[gene_id] = {"chrom": chrom,
                                  "start": int(start),
                                  "end": int(end),
                                  "strand": strand,
                                  "gene_biotype2": gene_biotype, 
                                  "gene_name": gene_name}
            else:
                if int(start) < genes[gene_id]["start"]:
                    genes[gene_id]["start"] = int(start)
                if int(end) > genes[gene_id]["end"]:
                    genes[gene_id]["end"] = int(end)

    if redo_genes:
        for gene in genes:
            d = genes[gene]
            chrom = d["chrom"]
            start = str(d["start"])
            end = str(d["end"])
            strand = d["strand"]
            gene_biotype2 = d["gene_biotype2"]
            gene_name = d["gene_name"]
            gene_biotype = assign_simplified_geneType(gene_biotype2)
            coordinates = ":".join([chrom, start, end, strand])
            gene_id = "|".join([gene, gene_name, gene_biotype, coordinates])

            gene_bed_line = "\t".join([chrom, start, end, gene_id, '0', strand]) + "\n"
            output_genes.write(gene_bed_line)

    if redo_exons or redo_genes:
        print("\r\tBEDs for {}... generated ".format(sp))


def bed_sort():
    files = dirs = os.listdir("./BEDs/")
    for file_name in files:
        print("\r\t"+file_name+"... sorting", end="")
        call = "sort -u -k1,1 -k2,2n -o '%s' '%s'"%("./BEDs/"+file_name,"./BEDs/"+file_name)
        subprocess.call(call, shell=True)
        print("\r\t"+file_name+"... sorted ")


def liftOver(arg_list):
	sp1 = arg_list[0]
	sp2 = arg_list[1]
	feature = arg_list[2]
	mM = arg_list[3]

	if not options.redo and (os.path.isfile("liftovers/{}to{}.{}.liftover".format(sp1, sp2, feature)) or os.path.isfile("liftovers/{}to{}.{}.unmapped".format(sp1, sp2, feature))):
		print("\t{} {} to {} already exists.. Skipping".format(sp1, feature, sp2))
	else:
		print("\t{} {} to {}... mapping".format(sp1, feature, sp2))
		oldFile = "BEDs/{}.{}.bed".format(sp1, feature)
		map_chain = [chain for chain in config_df.loc[config_df["assembly_version"]==sp1, "chainmap"].iloc[0].split(",") if sp2.lower() in chain.lower()][0]
		newFile = "liftovers/{}to{}.{}.liftover".format(sp1, sp2, feature)
		unMapped= "liftovers/{}to{}.{}.unmapped".format(sp1, sp2, feature)
		os.system("./scripts/liftOver -minMatch={} {} {} {} {}".format(mM/100, oldFile, map_chain, newFile, unMapped))
		print("\t{} {} to {}... done         ".format(sp1, feature, sp2))


def geneID_map_to_dict(file_name, dl = "\t"):
    dict_ = {}
    with open(file_name, 'r') as f:
        for line in f:
            line = line.strip().split(dl)
            dict_[line[0]] = {"gene_name": line[1],
                              "gene_type": line[2],
                              "coordinates": line[3],
                              "transcript_IDs": line[4].split("|")}
    return(dict_)


def parse_orthologs(line):
    geneA_ID = line[0]
    geneA_biotype = line[1]
    geneB_IDs = line[2].split(",")
    nexon = line[3].split(",")
    pcent = line[4].split(",")
    geneB_biotypes = line[5].split(",")

    return(geneA_ID, geneA_biotype, geneB_IDs, nexon, pcent, geneB_biotypes)

def count_biotypes(geneB_biotypes):
    tmp   = []
    none  = ["none"]
    pc    = ["protein_coding", "pc"]
    lncRNA = ["lncRNA", "TEC", "processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",
"non_coding", "sense_intronic", "sense_overlapping", "macro_lncRNA","antisense_RNA", 
"bidirectional_promoter_lncrna", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA"]

    ncRNA = ["ncRNA", "snRNA", "snoRNA", "rRNA", "rRNA_pseudogene", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "miRNA", "ribozyme", "sRNA", 
"scaRNA", "scRNA", "vault_RNA", "known_ncrna"]

    other = ["other", "processed_pseudogene", "transcribed_processed_pseudogene", "translated_processed_pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "transcribed_unitary_pseudogene", "polymorphic pseudogene", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "IG_C_pseudogene", "IG_D_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "IG_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "NOVEL"]

    for type_ in geneB_biotypes:
        if type_ in none:
            tmp.append("none")
        elif type_ in ncRNA:
            tmp.append("ncRNA")
        elif type_ in lncRNA:
            tmp.append("lncRNA")
        elif type_ in pc:
            tmp.append("pc")
        else:
            tmp.append("other")

    biotype_counts = [tmp.count("none"), tmp.count("ncRNA"), tmp.count("lncRNA"), tmp.count("pc"), tmp.count("other")]
    classification = classify_biotypes(biotype_counts)

    return(biotype_counts, classification)


def classify_biotypes(bcounts):
    classification = "other"
    #unique cases
    #  none              lncRNA            pc                other
    if bcounts[0]==1 and bcounts[1]==0 and bcounts[2]==0 and bcounts[3]==0:
        classification = "none"
    if bcounts[0]>=0 and bcounts[1]==1 and bcounts[2]==0 and bcounts[3]==0:
        classification = "lncRNA"
    if bcounts[0]>=0 and bcounts[1]==0 and bcounts[2]==1 and bcounts[3]==0:
        classification = "pc"
    if bcounts[0]>=0 and bcounts[1]==0 and bcounts[2]==0 and bcounts[3]==1:
        classification = "other"

    #multiple cases
    #  none              lncRNA            pc                other
    if bcounts[0]>1  and bcounts[1]==0 and bcounts[2]==0 and bcounts[3]==0:
        classification = "nones"
    if bcounts[0]>=0 and bcounts[1]>1  and bcounts[2]==0 and bcounts[3]==0:
        classification = "lncRNAs"
    if bcounts[0]>=0 and bcounts[1]==0 and bcounts[2]>1  and bcounts[3]==0:
        classification = "pcs"
    if bcounts[0]>=0 and bcounts[1]==0 and bcounts[2]==0 and bcounts[3]>1:
        classification = "others"

    #dual cases
    #  none              lncRNA            pc                other
    if bcounts[0]>=0 and bcounts[1]>=1 and bcounts[2]>=1 and bcounts[3]>=0:
        classification = "lncRNA_PC"
    if bcounts[0]>=0 and bcounts[1]>=1 and bcounts[2]==0 and bcounts[3]>=1:
        classification = "lncRNA_other"
    if bcounts[0]>=0 and bcounts[1]==0 and bcounts[2]>=1 and bcounts[3]>=1:
        classification = "PC_other"

    return(classification)


### Read and check config file
print("Checking Config_file...")

config_df = read_config("./config")
config_df.set_index("species", inplace = True)
config_df["default"] = "False|False"
for i in config_df.index:
    #Check if assembly_version exists in dictionaries
    if not config_df.loc[i]["assembly_version"].lower() in dictionaries["chain_maps"]:
        print("'%s' is not a valid assembly_version. Please use one of the following values for default analysis: %s"%(config_df.loc[i]["assembly_version"], ",".join(dictionaries["chain_maps"])))
        sys.exit(1)
    #Add and download default annotations in config if not provided
    if not config_df.loc[i]["annotation"]:
        if "danrer" in config_df.loc[i]["assembly_version"].lower():
        	default_gtf = dictionaries["gtfs_ensembl_r98"][config_df.loc[i]["assembly_version"].lower()].split("/")[-1]
        	download_default_files(dictionaries["gtfs_ensembl_r98"][config_df.loc[i]["assembly_version"].lower()], "GTFs")
        else:
        	default_gtf = dictionaries["gtfs"][config_df.loc[i]["assembly_version"].lower()].split("/")[-1]
        	download_default_files(dictionaries["gtfs"][config_df.loc[i]["assembly_version"].lower()], "GTFs")
        config_df.at[i, 'default'] = change_default("gtf", "True", config_df, i)
        config_df.at[i, 'annotation'] = "GTFs/"+default_gtf
    else:
        check_file([config_df.loc[i]["annotation"]])
    #Add and download default chainmaps in config if not provided
    if not config_df.loc[i]["chainmap"]:
        chainmaps = []
        config_df.at[i, 'default'] = change_default("chainmap", "True", config_df, i)
        for j in config_df.index:
            if i!=j:
                default_chainmap_path = dictionaries["chain_maps"][config_df.loc[i]["assembly_version"].lower()][config_df.loc[j]["assembly_version"].lower()]
                default_chainmap_name = "chainmaps/"+default_chainmap_path.split("/")[-1]
                chainmaps.append(default_chainmap_name)
                download_default_files(default_chainmap_path, "chainmaps")
        config_df.at[i, 'chainmap'] = ",".join(chainmaps)
    else:
        check_file(config_df.loc[i]["chainmap"].split(","))

print("Config_file is correct..")

### Generate GeneID - geneName - geneBiotype - transcripts maps
#check_folder("maps")
output_file = "maps/"
if os.path.exists(output_file):
    # If folder exists, delete all contents
    shutil.rmtree(output_file)
os.makedirs(output_file)
print("Generating geneID - geneName - geneBiotype - transcripts maps...")

gtfFiles_assembly = list(zip(list(config_df["annotation"]), list(config_df["assembly_version"])))
pool = mp.Pool(mp.cpu_count())
pool.map(generate_maps, gtfFiles_assembly)
pool.close()


### Generate BED with genes/exons from GTF files
#check_folder("BEDs")
output_file = "BEDs/"
if os.path.exists(output_file):
    # If folder exists, delete all contents
    shutil.rmtree(output_file)
os.makedirs(output_file)
print("Generating BED files for exons and genes...")

gtfFiles_assembly_level = list( zip( list(config_df["annotation"]), list(config_df["assembly_version"]), [gene_level]*len(list(config_df["assembly_version"])) ) )
pool = mp.Pool(mp.cpu_count())
pool.map(generate_beds, gtfFiles_assembly_level)
pool.close()

# Remove possible exon duplications due to different isoforms in GTF files
for fileName in os.listdir("BEDs"):
    call = 'sort %s | uniq > %s; mv %s %s'%('BEDs/'+fileName, 'BEDs/tmp.bed', 'BEDs/tmp.bed', 'BEDs/'+fileName)
    subprocess.call(call, shell=True, executable='/bin/bash')

### LifOver exons/genes

#check_folder("liftovers/")
output_file = "liftovers/"
if os.path.exists(output_file):
    # If folder exists, delete all contents
    shutil.rmtree(output_file)
os.makedirs(output_file)

print("LiftOver features between species...")

features = ["exons", "genes"] if gene_level else ["exons"]
arguments = []
for sp1 in config_df["assembly_version"]:
    for sp2 in config_df["assembly_version"]:
        if sp1 == sp2: continue
        for feature in features:
            arg_list = [sp1, sp2, feature, minMatch]
            arguments.append(arg_list)

pool = mp.Pool(mp.cpu_count())
pool.map(liftOver, arguments)
pool.close()

### Intersect LiftOvers
#check_folder("overlaps/")
output_file = "overlaps/"
if os.path.exists(output_file):
    # If folder exists, delete all contents
    shutil.rmtree(output_file)
os.makedirs(output_file)

print("Intersecting LiftOver...")
for i in config_df.index:
	for j in config_df.index:
		if i == j: continue
		for f in ["exons", "genes"]:
			if f == "genes" and not gene_level: continue
			sp1=config_df.loc[i]["assembly_version"].lower()
			sp2=config_df.loc[j]["assembly_version"].lower()

			# Overlap exons_span with gene_span
			if os.path.isfile('overlaps/%sto%s.%s.overlap'%(sp1, sp2, 'gene')):
				pass
			else:
				#output_exon = 'overlaps/%sto%s.%s.overlap'%(sp1, sp2, 'exon')
				output_exonStrand = 'overlaps/%sto%s.%s.overlap'%(sp1, sp2, 'exonStrand')
				#output_gene = 'overlaps/%sto%s.%s.overlap'%(sp1, sp2, 'gene')
				output_geneStrand = 'overlaps/%sto%s.%s.overlap'%(sp1, sp2, 'geneStrand')

				lifover_input = 'liftovers/%sto%s.%s.liftover'%(sp1, sp2, 'exons')
				bed_exonsInput = 'BEDs/%s.%s.bed'%(sp2, 'exons')
				bed_genesInput = 'BEDs/%s.%s.bed'%(sp2, 'genes')

				#call_exon = 'intersectBed -wao -a %s -b %s > %s'%(lifover_input, bed_exonsInput, output_exon)
				call_exonStrand = 'intersectBed -wao -s -a %s -b %s > %s'%(lifover_input, bed_exonsInput, output_exonStrand)
				#call_gene = 'intersectBed -wao -a %s -b %s > %s'%(lifover_input, bed_genesInput, output_gene)
				call_geneStrand = 'intersectBed -wao -s -a %s -b %s > %s'%(lifover_input, bed_genesInput, output_geneStrand)

				for c in [call_exonStrand, call_geneStrand]:
					subprocess.call(c, shell=True, executable='/bin/bash')
				#for c in [call_exon, call_exonStrand, call_gene, call_geneStrand]:
				#	subprocess.call(c, shell=True, executable='/bin/bash')

			# Overlap exons with exons and genes with genes
			output = 'overlaps/%sto%s.%s.overlap'%(sp1, sp2, f)
			if not options.redo and os.path.isfile(output):
				print("\t{} {} to {} already exist.. Skipping".format(i, f, j))
			else:
				print("\r\t{} {} to {}... intersecting".format(i, f, j), end="")
				lifover_input = 'liftovers/%sto%s.%s.liftover'%(sp1, sp2, f)
				bed_input = 'BEDs/%s.%s.bed'%(sp2, f)
				call = 'intersectBed -wao -s -a %s -b %s > %s'%(lifover_input, bed_input, output)
				subprocess.call(call, shell=True, executable='/bin/bash')
				print("\r\t{} {} to {}... done         ".format(i, f, j))


### Parse Overlap
#check_folder("overlaps/")
output_file = "orthology/"
if os.path.exists(output_file):
    # If folder exists, delete all contents
    shutil.rmtree(output_file)
os.makedirs(output_file)
print("Predicting orthologues...")

for i in config_df.index:
	for j in config_df.index:

		if i == j: continue

		sp1=config_df.loc[i]["assembly_version"].lower()
		sp2=config_df.loc[j]["assembly_version"].lower()

		sp1_gID_Name_Type = geneID_map_to_dict('maps/%s.geneID_map.txt'%(sp1))
		sp2_gID_Name_Type = geneID_map_to_dict('maps/%s.geneID_map.txt'%(sp2))

		for f in ["exonStrand", "genes"]:
			if f == "genes" and not gene_level: continue
			if not options.redo and os.path.isfile('orthology/%sto%s.%s'%(sp1, sp2, f)):
				print("\t{} {} to {} orthology already exist.. Skipping".format(i, f, j))
				continue
			print("\r\t{} {} to {}... finding orthologues".format(i, f, j), end="")
			overlaps = "overlaps/%sto%s.%s.overlap"%(sp1, sp2, f)

			genes = {}
			for line in open(overlaps, 'r'):
				# read line overlap between sp1 and sp2
				line = line.strip().split("\t")
				exonA = line[3]
				geneA_ID = exonA.split("|")[0]
				exonB = line[9]
				geneB_ID = exonB.split("|")[0]

				# calculate overlapping %
				lA = int(line[2])-int(line[1]) #length geneA
				lB = int(line[8])-int(line[7]) #length geneB
				o = int(line[12])              #overlap

				rA = (o*100)/lA                    #overlap of A to B
				rB = 0 if lB == 0 else (o*100)/lB  #overlap of B to A

				# Harmonic mean
				try:
					hm = 2.0/((1.0/rA)+(1.0/rB))
				except ZeroDivisionError:
					hm = 0

				# Gene_Name and Gene_Type
				geneAtype = "none" if geneA_ID == "." else sp1_gID_Name_Type[geneA_ID]["gene_type"]
				geneBtype = "none" if geneB_ID == "." else sp2_gID_Name_Type[geneB_ID]["gene_type"]

				# Add sp1 gene overlaps to dictionary
				if not geneA_ID in genes:
					genes[geneA_ID] = {"geneAtype": geneAtype,
									   "overlaps": [],
									   "geneBtypes": [],
									   "exons": {}}

				if not geneB_ID in genes[geneA_ID]["overlaps"]:
					genes[geneA_ID]["overlaps"].append(geneB_ID)
					genes[geneA_ID]["geneBtypes"].append(geneBtype)

				# Add sp1 exon overlaps to dictionary
				if not exonA in genes[geneA_ID]["exons"]:
					genes[geneA_ID]["exons"][exonA] = {}

				if not geneB_ID in genes[geneA_ID]["exons"][exonA]:
					genes[geneA_ID]["exons"][exonA][geneB_ID] = [rA, rB, hm]

				# Update exonB if new exon have highr A>B overlap
				# NOTE ---> and higher B>A overlap? change % of overlap with HM ??
				if rA > genes[geneA_ID]["exons"][exonA][geneB_ID][0]:
					genes[geneA_ID]["exons"][exonA][geneB_ID] = [rA, rB, hm]

			# Summarize overlaps into output
			output_file = open('orthology/%sto%s.%s'%(sp1, sp2, f), "w")
			header = "\t".join(["GeneA_ID", "GeneA_biotype", "GeneB_IDs", "GeneB_biotypes", "GeneB_Classification", "GeneA_number_overlapping_exons", "GeneA_percent_overlapping"])
			output_file.write(header+"\n")
			for geneA_ID in genes:
				overlaps = genes[geneA_ID]["overlaps"]
				geneAtype = genes[geneA_ID]["geneAtype"]
				geneBtype = genes[geneA_ID]["geneBtypes"]
				number_of_exons = []
				percentageA = []
				percentageB = []
				for geneB_ID in overlaps:
					n = 0
					mA, mB = 0, 0
					for exon in genes[geneA_ID]["exons"]:
						if geneB_ID in genes[geneA_ID]["exons"][exon]:
							n+=1
							# %overlap A>B (rA)
							mA+=genes[geneA_ID]["exons"][exon][geneB_ID][0]
							mB+=genes[geneA_ID]["exons"][exon][geneB_ID][1]
					number_of_exons.append(str(n))
					percentageA.append(str(mA/n))
					percentageB.append(str(mB/n))

				percentageA = [ float(x) for x in percentageA ]
				percentageB = [ float(x) for x in percentageB ]
				number_of_exons = [ int(x) for x in number_of_exons ]

				percentageA, number_of_exons, overlaps, geneBtype = zip(*sorted(zip(percentageA, number_of_exons, overlaps, geneBtype), reverse=True))
				if "." in overlaps and len(overlaps)>1:
					overlaps = list(overlaps)
					geneBtype = list(geneBtype)
					overlaps.remove('.')
					geneBtype.remove("none")
					overlaps = tuple(overlaps)
					geneBtype = tuple(geneBtype)
				percentageA = [ str(x) for x in percentageA ]
				number_of_exons = [ str(x) for x in number_of_exons ]
				biotype_counts, classification = count_biotypes(geneBtype)
				#biotype_counts, classification = count_biotypes(geneAtype)
				new_line = "\t".join([geneA_ID, geneAtype, ",".join(overlaps), ",".join(geneBtype), classification, ",".join(number_of_exons), ",".join(percentageA)])
				output_file.write(new_line+"\n")

			output_file.close()
			print("\r\t{} {} to {}... done               ".format(i, f, j))


###########################################################################################################################
###########################################################################################################################
###############################                                                        ####################################
###############################               Bi-Directional Start                     ####################################                  
###############################                                                        ####################################
###########################################################################################################################
###########################################################################################################################


input_path = './orthology'
files = [x for x in os.listdir(input_path)]

one_to_one_sp1_A_B = set()
one_to_one_sp2_B_A = set()
all_sp1_sp2 = set()


# Process each input file
for input_file2 in files:
    if gene_level == True:
        if input_file2.endswith(".exonStrand"):
            continue
        cleaned_filename_ortho = input_file2.replace('.genes', '')  # Remove '.genes'
        values1 = cleaned_filename_ortho.strip().split('to')  # Split at 'to'
        sp_name = values1[0]
        input_file = input_file2

    if gene_level == False:
        if input_file2.endswith(".genes"):
            continue
        cleaned_filename_ortho = input_file2.replace('.exonStrand', '')  # Remove '.exonStrand'
        values1 = cleaned_filename_ortho.strip().split('to')  # Split at 'to'
        sp_name = values1[0]
        input_file = input_file2
    #print(input_path, input_file, sp_name, cleaned_filename_ortho)
    with open(os.path.join(input_path, input_file), 'r') as handle:
        next(handle)  # Skip header line

        for line in handle:
            line = line.strip().split("\t")
            gene_id = line[0]
            gene_nameA = line[1]
            geneB_id_elements = line[2].split(",")
            geneB_name_elements = line[3].split(",")
            gene_info_A = ''
            gene_info_B = ''
            for geneB_id , geneB_name in zip(geneB_id_elements, geneB_name_elements):

                if sp1 == sp_name:
                    gene_info_A = (gene_id, gene_nameA, geneB_id, geneB_name)
                    one_to_one_sp1_A_B.add(gene_info_A)
                    all_sp1_sp2.add(gene_info_A)
                    if geneB_id != '.':
                        gene_info_B = (geneB_id, geneB_name, gene_id, gene_nameA)
                        all_sp1_sp2.add(gene_info_B)


                if sp2 == sp_name:
                    gene_info_B = (gene_id, gene_nameA, geneB_id, geneB_name)
                    one_to_one_sp2_B_A.add(gene_info_B)
                    all_sp1_sp2.add(gene_info_B)
                    if geneB_id != '.':
                        gene_info_A = (geneB_id, geneB_name, gene_id, gene_nameA)
                        all_sp1_sp2.add(gene_info_A)

# Creating mappings
one_to_one_sp1_A_B_map = {}
for gene_id_A1, gene_biotype_A1, gene_id_B1, gene_biotype_B1 in one_to_one_sp1_A_B:
    if gene_id_A1 not in one_to_one_sp1_A_B_map:
        one_to_one_sp1_A_B_map[gene_id_A1] = []
    one_to_one_sp1_A_B_map[gene_id_A1].append((gene_id_A1, gene_biotype_A1, gene_id_B1, gene_biotype_B1))

one_to_one_sp2_B_A_map = {}
for gene_id_A2, gene_biotype_A2, gene_id_B2, gene_biotype_B2 in one_to_one_sp2_B_A:
    if gene_id_A2 not in one_to_one_sp2_B_A_map:
        one_to_one_sp2_B_A_map[gene_id_A2] = []
    one_to_one_sp2_B_A_map[gene_id_A2].append((gene_id_A2, gene_biotype_A2, gene_id_B2, gene_biotype_B2))

all_sp1_sp2_map = {}
for gene_id_A3, gene_biotype_A3, gene_id_B3, gene_biotype_B3 in all_sp1_sp2:
    if gene_id_A3 not in all_sp1_sp2_map:
        all_sp1_sp2_map[gene_id_A3] = []
    all_sp1_sp2_map[gene_id_A3].append((gene_id_A3, gene_biotype_A3, gene_id_B3, gene_biotype_B3))


############################################################################################################

statistic_2l = "statistics/lncRNA/"
if os.path.exists(statistic_2l):
    shutil.rmtree(statistic_2l)
os.makedirs(statistic_2l)

statistic_2p = "statistics/protein/"
if os.path.exists(statistic_2p):
    shutil.rmtree(statistic_2p)
os.makedirs(statistic_2p)

Orthologous_l2 = "Orthologous/lncRNA/"
if os.path.exists(Orthologous_l2):
    shutil.rmtree(Orthologous_l2)
os.makedirs(Orthologous_l2)

Orthologous_p2 = "Orthologous/protein/"
if os.path.exists(Orthologous_p2):
    shutil.rmtree(Orthologous_p2)
os.makedirs(Orthologous_p2)

Bi_Directional_plot = "plots/"
if os.path.exists(Bi_Directional_plot):
    shutil.rmtree(Bi_Directional_plot)
os.makedirs(Bi_Directional_plot)

input_path = "BEDs/"
files = os.listdir(input_path)
#######################################################################
for input_file2 in files:
    if gene_level == True:
        if input_file2.endswith("exons.bed"):
            continue
        cleaned_filename = input_file2.replace('.genes.bed', '')  # Remove '.genes'
        values1 = cleaned_filename.strip().split('to')  # Split at 'to'
        sp_name = values1[0]
        input_file = input_file2
        file_path = os.path.join(input_path, input_file)


    if gene_level == False:
        if input_file2.endswith(".genes.bed"):
            continue
        cleaned_filename = input_file2.replace('.exons.bed', '')  # Remove '.exonStrand'
        values1 = cleaned_filename.strip().split('to')  # Split at 'to'
        sp_name = values1[0]
        input_file = input_file2
        file_path = os.path.join(input_path, input_file)
    #print(input_path, input_file, sp_name, cleaned_filename)
    with open(file_path, 'r') as handle:
        
        # Initialize data structures
        total_pc = set()
        not_lifted_pc = set()
        lifted_pc = set()
        lifted_one_none_pc = set()
        lifted_one_one_pc = set()
        lifted_one_many_pc = set()
        lifted_many_many_pc = set()
        lifted_one_one_pc_lncRNA = set()
        lifted_one_one_pc_pc = set()
        lifted_one_one_pc_other = set()
        lifted_one_half_pc = set()
        ortholog_one_many_pc = set()
        ortholog_many_many_pc = set()

        total_lncRNA = set()
        not_lifted_lncRNA = set()
        lifted_lncRNA = set()
        lifted_one_none_lncRNA = set()
        lifted_one_one_lncRNA = set()
        lifted_one_half_lncRNA = set()
        lifted_one_many_lncRNA = set()
        lifted_many_many_lncRNA = set()
        lifted_one_one_lncRNA_lncRNA = set()
        lifted_one_one_lncRNA_pc = set()
        lifted_one_one_lncRNA_other = set()
        ortholog_one_many_lncRNA = set()
        ortholog_many_many_lncRNA = set()

        data = set()

        # Parse the input file
        for line in handle:
            line = line.strip().split("\t")
            exonID1 = line[3]
            gene_id = exonID1.split("|")[0]
            gene_biotype = exonID1.split("|")[2]
            data.add((gene_id, gene_biotype))
            
        
        for gene_id, gene_biotype in data:
            if gene_biotype == 'pc':
                total_pc.add(gene_id)
            elif gene_biotype == 'lncRNA':
                total_lncRNA.add(gene_id)

            gene_id_A_half_no = 0
            gene_id_B_half_no = 0
            gene_id_A_unic = set()
            gene_id_B_unic = set()
            gene_id_A_unic_no = 0
            gene_id_B_unic_no = 0

#######################################################################################################################################
#################                                                                                                   ###################
#################                                       SP1  start                                                  ###################
#################                                                                                                   ###################
#######################################################################################################################################
            #print(gene_id, gene_biotype)
            if gene_id in all_sp1_sp2_map and sp1 == cleaned_filename:
                cleaned_filename2= sp1 + "to" + sp2
                for gene_id_A, gene_biotype_A, gene_id_B, gene_biotype_B in all_sp1_sp2_map[gene_id]:
                    gene_id_A_unic.add((gene_id_A, gene_biotype_A))

                    if gene_id_B != '.':
                        gene_id_B_unic.add((gene_id_B, gene_biotype_B))

                    if gene_id_A in one_to_one_sp1_A_B_map:
                        for gene_id1_A, gene_biotype1_A, gene_id1_B, gene_biotype1_B in one_to_one_sp1_A_B_map[gene_id_A]:
                            gene_id_A_unic.add((gene_id1_A, gene_biotype1_A))
                            if gene_id1_B != '.':
                                gene_id_A_half_no = 1
                                gene_id_B_unic.add((gene_id1_B, gene_biotype1_B))


                    if gene_id_B in one_to_one_sp2_B_A_map:
                        for gene_id2_B, gene_biotype2_B, gene_id2_A, gene_biotype2_A in one_to_one_sp2_B_A_map[gene_id_B]:
                            gene_id_B_unic.add((gene_id2_B, gene_biotype2_B))
                            if gene_id_A != '.':
                                gene_id_B_half_no = 1
                                gene_id_A_unic.add((gene_id2_A, gene_biotype2_A))

                    if gene_id_B in all_sp1_sp2_map:
                        for gene_id3_B, gene_biotype3_B, gene_id3_A, gene_biotype3_A in all_sp1_sp2_map[gene_id_B]:
                            gene_id_B_unic.add((gene_id3_B, gene_biotype3_B))
                            if gene_id3_A != '.':
                                gene_id_A_unic.add((gene_id3_A, gene_biotype3_A))
         

                gene_id_A_unic_no = len(gene_id_A_unic)
                gene_id_B_unic_no = len(gene_id_B_unic)

#################################################################################################################
#################################################################################################################

                ################################################################################
                # One-to-one orthology
                ################################################################################

                if gene_id_A_unic_no == 1 and gene_id_B_unic_no == 1 and gene_id_A_half_no == 1 and gene_id_B_half_no == 1:
                    for (gene_idA, gene_biotypeA), (gene_idB, gene_biotypeB) in zip(gene_id_A_unic, gene_id_B_unic):

                        if gene_biotypeA == 'lncRNA':
                            lifted_one_one_lncRNA.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))

                            if gene_biotypeB == 'lncRNA':
                                lifted_one_one_lncRNA_lncRNA.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_lncRNA.add(gene_id)

                            elif gene_biotypeB == 'pc':
                                lifted_one_one_lncRNA_pc.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_lncRNA.add(gene_id)

                            elif gene_biotypeB in ['ncRNA', 'other']:
                                lifted_one_one_lncRNA_other.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_lncRNA.add(gene_id)

                        if gene_biotypeA == 'pc':
                            lifted_one_one_pc.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))

                            if gene_biotypeB == 'pc':
                                lifted_one_one_pc_pc.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_pc.add(gene_id)

                            elif gene_biotypeB == 'lncRNA':
                                lifted_one_one_pc_lncRNA.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_pc.add(gene_id)

                            elif gene_biotypeB in ['ncRNA', 'other']:
                                lifted_one_one_pc_other.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_pc.add(gene_id)

                ################################################################################
                # One-to-half orthology
                ################################################################################
                elif gene_id_A_unic_no == 1 and gene_id_B_unic_no == 1:

                    if gene_id_A_half_no == 1 and gene_id_B_half_no == 0:
                        for (gene_idA, gene_biotypeA), (gene_idB, gene_biotypeB) in zip(gene_id_A_unic, gene_id_B_unic):
                            if gene_biotypeA == 'lncRNA':
                                lifted_one_half_lncRNA.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_lncRNA.add(gene_id)

                            if gene_biotypeA == 'pc':
                                lifted_one_half_pc.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_pc.add(gene_id)


                ################################################################################
                # One-to-none orthology
                ################################################################################
                elif gene_id_A_unic_no == 1 and gene_id_B_unic_no == 0:
                    for gene_idA, gene_biotypeA in gene_id_A_unic:
                        if gene_biotype == 'lncRNA':
                            lifted_one_none_lncRNA.add(gene_idA)
                            lifted_lncRNA.add(gene_idA)


                        if gene_biotype == 'pc':
                            lifted_one_none_pc.add(gene_idA)
                            lifted_pc.add(gene_idA)
                ################################################################################
                # One-to-many orthology
                ################################################################################
                elif gene_id_A_unic_no == 1 and gene_id_B_unic_no >= 2:
                    if gene_biotype == 'lncRNA':
                        gene_id_A_str = ' '.join([f"{id_} {biotype}" for id_, biotype in gene_id_A_unic])
                        ortholog_one_many_lncRNA.add(f"{gene_id_A_str} {tuple(gene_id_B_unic)}")
                        lifted_one_many_lncRNA.add(gene_id)
                        lifted_lncRNA.add(gene_id)

                    if gene_biotype == 'pc':
                        gene_id_A_str = ' '.join([f"{id_} {biotype}" for id_, biotype in gene_id_A_unic])
                        ortholog_one_many_pc.add(f"{gene_id_A_str} {tuple(gene_id_B_unic)}")
                        lifted_one_many_pc.add(gene_id)
                        lifted_pc.add(gene_id)

                elif gene_id_A_unic_no >= 2 and gene_id_B_unic_no == 1:

                    if gene_biotype == 'lncRNA':
                        gene_id_A_str = ' '.join([f"{id_} {biotype}" for id_, biotype in gene_id_A_unic])
                        ortholog_one_many_lncRNA.add(f"{gene_id_A_str} {tuple(gene_id_B_unic)}")
                        lifted_one_many_lncRNA.add(gene_id)
                        lifted_lncRNA.add(gene_id)

                    if gene_biotype == 'pc':
                        gene_id_A_str = ' '.join([f"{id_} {biotype}" for id_, biotype in gene_id_A_unic])
                        ortholog_one_many_pc.add(f"{gene_id_A_str} {tuple(gene_id_B_unic)}")
                        lifted_one_many_pc.add(gene_id)
                        lifted_pc.add(gene_id)
                    ################################################################################
                    # Many-to-many orthology
                    ################################################################################
                elif gene_id_A_unic_no >= 2 and gene_id_B_unic_no >= 2:

                    if gene_biotype == 'lncRNA':
                        ortholog_many_many_lncRNA.add(f"{tuple(gene_id_A_unic)} {tuple(gene_id_B_unic)}")
                        lifted_many_many_lncRNA.add(gene_id)
                        lifted_lncRNA.add(gene_id)

                    if gene_biotype == 'pc':
                        ortholog_many_many_pc.add(f"{tuple(gene_id_A_unic)} {tuple(gene_id_B_unic)}")
                        lifted_many_many_pc.add(gene_id)
                        lifted_pc.add(gene_id)

#######################################################################################################################################
#################                                                                                                   ###################
#################                                       SP1  end                                                  ###################
#################                                                                                                   ###################
#######################################################################################################################################

#######################################################################################################################################
#################                                                                                                   ###################
#################                                       SP2  start                                                  ###################
#################                                                                                                   ###################
#######################################################################################################################################

            if gene_id in all_sp1_sp2_map and sp2 == cleaned_filename:
                cleaned_filename2= sp2 + "to" + sp1
                for gene_id_A, gene_biotype_A, gene_id_B, gene_biotype_B in all_sp1_sp2_map[gene_id]:
                    gene_id_A_unic.add((gene_id_A, gene_biotype_A))

                    if gene_id_B != '.':
                        gene_id_B_unic.add((gene_id_B, gene_biotype_B))

                    if gene_id_A in one_to_one_sp2_B_A_map:
                        for gene_id2_A, gene_biotype2_A, gene_id2_B, gene_biotype2_B in one_to_one_sp2_B_A_map[gene_id_A]:
                            gene_id_A_unic.add((gene_id2_A, gene_biotype2_A))
                            if gene_id2_B != '.':
                                gene_id_A_half_no = 1
                                gene_id_B_unic.add((gene_id2_B, gene_biotype2_B))

                    if gene_id_B in one_to_one_sp1_A_B_map:
                        for gene_id1_B, gene_biotype1_B, gene_id1_A, gene_biotype1_A in one_to_one_sp1_A_B_map[gene_id_B]:
                            gene_id_B_unic.add((gene_id1_B, gene_biotype1_B))
                            if gene_id1_A != '.':
                                gene_id_B_half_no = 1
                                gene_id_A_unic.add((gene_id1_A, gene_biotype1_A))


                    if gene_id_B in all_sp1_sp2_map:
                        for gene_id3_B, gene_biotype3_B, gene_id3_A, gene_biotype3_A in all_sp1_sp2_map[gene_id_B]:
                            gene_id_B_unic.add((gene_id3_B, gene_biotype3_B))
                            if gene_id3_A != '.':
                                gene_id_A_unic.add((gene_id3_A, gene_biotype3_A))


                gene_id_A_unic_no = len(gene_id_A_unic)
                gene_id_B_unic_no = len(gene_id_B_unic)

                ################################################################################
                # One-to-one orthology
                ################################################################################
                if gene_id_A_unic_no == 1 and gene_id_B_unic_no == 1 and gene_id_A_half_no == 1 and gene_id_B_half_no == 1:
                    for (gene_idA, gene_biotypeA), (gene_idB, gene_biotypeB) in zip(gene_id_A_unic, gene_id_B_unic):
                        if gene_biotypeA == 'lncRNA':
                            lifted_one_one_lncRNA.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))

                            if gene_biotypeB == 'lncRNA':
                                lifted_one_one_lncRNA_lncRNA.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_lncRNA.add(gene_id)

                            elif gene_biotypeB == 'pc':
                                lifted_one_one_lncRNA_pc.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_lncRNA.add(gene_id)

                            elif gene_biotypeB in ['ncRNA', 'other']:
                                lifted_one_one_lncRNA_other.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_lncRNA.add(gene_id)


                        if gene_biotypeA == 'pc':
                            lifted_one_one_pc.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))

                            if gene_biotypeB == 'pc':
                                lifted_one_one_pc_pc.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_pc.add(gene_id)

                            elif gene_biotypeB == 'lncRNA':
                                lifted_one_one_pc_lncRNA.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_pc.add(gene_id)

                            elif gene_biotypeB in ['ncRNA', 'other']:
                                lifted_one_one_pc_other.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_pc.add(gene_id)

                ################################################################################
                # One-to-half orthology
                ################################################################################
                elif gene_id_A_unic_no == 1 and gene_id_B_unic_no == 1:

                    if gene_id_A_half_no == 1 and gene_id_B_half_no == 0:
                        for (gene_idA, gene_biotypeA), (gene_idB, gene_biotypeB) in zip(gene_id_A_unic, gene_id_B_unic):
                            if gene_biotypeA == 'lncRNA':
                                lifted_one_half_lncRNA.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_lncRNA.add(gene_id)

                            if gene_biotypeA == 'pc':
                                lifted_one_half_pc.add((gene_idA, gene_biotypeA, gene_idB, gene_biotypeB))
                                lifted_pc.add(gene_id)


                ################################################################################
                # One-to-none orthology
                ################################################################################
                elif gene_id_A_unic_no == 1 and gene_id_B_unic_no == 0:
                    for gene_idA, gene_biotypeA in gene_id_A_unic:
                        if gene_biotype == 'lncRNA':
                            lifted_one_none_lncRNA.add(gene_idA)
                            lifted_lncRNA.add(gene_idA)

                        if gene_biotype == 'pc':
                            lifted_one_none_pc.add(gene_idA)
                            lifted_pc.add(gene_idA)

                ################################################################################
                # One-to-many orthology
                ################################################################################
                elif gene_id_A_unic_no == 1 and gene_id_B_unic_no >= 2:
                    if gene_biotype == 'lncRNA':
                        gene_id_A_str = ' '.join([f"{id_} {biotype}" for id_, biotype in gene_id_A_unic])
                        ortholog_one_many_lncRNA.add(f"{gene_id_A_str} {tuple(gene_id_B_unic)}")
                        lifted_one_many_lncRNA.add(gene_id)
                        lifted_lncRNA.add(gene_id)

                    if gene_biotype == 'pc':
                        gene_id_A_str = ' '.join([f"{id_} {biotype}" for id_, biotype in gene_id_A_unic])
                        ortholog_one_many_pc.add(f"{gene_id_A_str} {tuple(gene_id_B_unic)}")
                        lifted_one_many_pc.add(gene_id)
                        lifted_pc.add(gene_id)

                elif gene_id_A_unic_no >= 2 and gene_id_B_unic_no == 1:
                    if gene_biotype == 'lncRNA':
                        gene_id_A_str = ' '.join([f"{id_} {biotype}" for id_, biotype in gene_id_A_unic])
                        ortholog_one_many_lncRNA.add(f"{gene_id_A_str} {tuple(gene_id_B_unic)}")
                        lifted_one_many_lncRNA.add(gene_id)
                        lifted_lncRNA.add(gene_id)

                    if gene_biotype == 'pc':
                        gene_id_A_str = ' '.join([f"{id_} {biotype}" for id_, biotype in gene_id_A_unic])
                        ortholog_one_many_pc.add(f"{gene_id_A_str} {tuple(gene_id_B_unic)}")
                        lifted_one_many_pc.add(gene_id)
                        lifted_pc.add(gene_id)

                ################################################################################
                # Many-to-many orthology
                ################################################################################
                elif gene_id_A_unic_no >= 2 and gene_id_B_unic_no >= 2:
                    if gene_biotype == 'lncRNA':
                        ortholog_many_many_lncRNA.add(f"{tuple(gene_id_A_unic)} {tuple(gene_id_B_unic)}")
                        lifted_many_many_lncRNA.add(gene_id)
                        lifted_lncRNA.add(gene_id)

                    if gene_biotype == 'pc':
                        ortholog_many_many_pc.add(f"{tuple(gene_id_A_unic)} {tuple(gene_id_B_unic)}")
                        lifted_many_many_pc.add(gene_id)
                        lifted_pc.add(gene_id)

#######################################################################################################################################
#################                                                                                                   ###################
#################                                       SP2  end                                                  ###################
#################                                                                                                   ###################
#######################################################################################################################################

        total_lncRNA_no = len(total_lncRNA)

        lifted_lncRNA_no = len(lifted_lncRNA)

        not_lifted_lncRNA_no = total_lncRNA_no - lifted_lncRNA_no

        lifted_one_half_lncRNA_no = len(lifted_one_half_lncRNA)

        lifted_one_many_lncRNA_no = len(lifted_one_many_lncRNA)

        lifted_many_many_lncRNA_no = len(lifted_many_many_lncRNA)

        lifted_one_none_lncRNA_no = len(lifted_one_none_lncRNA)

        lifted_one_one_lncRNA_no = len(lifted_one_one_lncRNA)

        lifted_one_one_lncRNA_lncRNA_no = len(lifted_one_one_lncRNA_lncRNA)

        lifted_one_one_lncRNA_pc_no = len(lifted_one_one_lncRNA_pc)

        lifted_one_one_lncRNA_other_no = len(lifted_one_one_lncRNA_other)

        total_pc_no = len(total_pc)

        lifted_pc_no = len(lifted_pc)

        not_lifted_pc_no = total_pc_no - lifted_pc_no

        lifted_one_half_pc_no = len(lifted_one_half_pc)

        lifted_one_many_pc_no = len(lifted_one_many_pc)

        lifted_many_many_pc_no = len(lifted_many_many_pc)

        lifted_one_none_pc_no = len(lifted_one_none_pc)

        lifted_one_one_pc_no = len(lifted_one_one_pc)

        lifted_one_one_pc_pc_no = len(lifted_one_one_pc_pc)

        lifted_one_one_pc_lncRNA_no = len(lifted_one_one_pc_lncRNA)

        lifted_one_one_pc_other_no = len(lifted_one_one_pc_other)

###########################################################################################
###########################################################################################

        one_lncRNA_to_one_lncRNA_out2 = os.path.join(Orthologous_l2, "one_lncRNA_to_one_lncRNA_" + cleaned_filename2 + ".txt")
        one_lncRNA_half_out2 = os.path.join(Orthologous_l2, "one_lncRNA_half_" + cleaned_filename2 + ".txt")
        lncRNA_one_to_one_out2 = os.path.join(Orthologous_l2, "one_to_one_" + cleaned_filename2 + ".txt")
        lncRNA_one_to_many_out2 = os.path.join(Orthologous_l2, "one_to_many_" + cleaned_filename2 + ".txt")
        lncRNA_many_to_many_out2 = os.path.join(Orthologous_l2, "many_to_many_" + cleaned_filename2 + ".txt")


        with open(one_lncRNA_to_one_lncRNA_out2, 'w') as output_file:
            for gene_info in lifted_one_one_lncRNA_lncRNA:
                output_file.write(f"{gene_info[0]}\t{gene_info[1]}\t{gene_info[2]}\t{gene_info[3]}\n")

        with open(lncRNA_one_to_one_out2, 'w') as output_file:
            for gene_info in lifted_one_one_lncRNA:
                output_file.write(f"{gene_info[0]}\t{gene_info[1]}\t{gene_info[2]}\t{gene_info[3]}\n")

        with open(one_lncRNA_half_out2, 'w') as output_file:
            for gene_info in lifted_one_half_lncRNA:
                output_file.write(f"{gene_info[0]}\t{gene_info[1]}\t{gene_info[2]}\t{gene_info[3]}\n")

        with open(lncRNA_one_to_many_out2, 'w') as output_file:
            for gene_info in ortholog_one_many_lncRNA:
                output_file.write(f"{gene_info}\n")

        with open(lncRNA_many_to_many_out2, 'w') as output_file:
            for gene_info in ortholog_many_many_lncRNA:
                output_file.write(f"{gene_info}\n")
#############################################################

        one_protein_to_one_protein_out2 = os.path.join(Orthologous_p2, "one_protein_to_one_protein_" + cleaned_filename2 + ".txt")
        pc_one_to_one_out2 = os.path.join(Orthologous_p2, "one_to_one_" + cleaned_filename2 + ".txt")
        pc_one_to_half_out2 = os.path.join(Orthologous_p2, "one_to_half_" + cleaned_filename2 + ".txt")
        pc_one_to_many_out2 = os.path.join(Orthologous_p2, "one_to_many_" + cleaned_filename2 + ".txt")
        pc_many_to_many_out2 = os.path.join(Orthologous_p2, "many_to_many_" + cleaned_filename2 + ".txt")

        with open(one_protein_to_one_protein_out2, 'w') as output_file:
            for gene_info in lifted_one_one_pc_pc:
                output_file.write(f"{gene_info[0]}\t{gene_info[1]}\t{gene_info[2]}\t{gene_info[3]}\n")

        with open(pc_one_to_one_out2, 'w') as output_file:
            for gene_info in lifted_one_one_pc:
                output_file.write(f"{gene_info[0]}\t{gene_info[1]}\t{gene_info[2]}\t{gene_info[3]}\n")

        with open(pc_one_to_half_out2, 'w') as output_file:
            for gene_info in lifted_one_half_pc:
                output_file.write(f"{gene_info[0]}\t{gene_info[1]}\t{gene_info[2]}\t{gene_info[3]}\n")

        with open(pc_one_to_many_out2, 'w') as output_file:
            for gene_info in ortholog_one_many_pc:
                output_file.write(f"{gene_info}\n")

        with open(pc_many_to_many_out2, 'w') as output_file:
            for gene_info in ortholog_many_many_pc:
                output_file.write(f"{gene_info}\n")
###########################################################################################


        no_lncRNA2 = 0
        lifted_spesi_data_lncRNA = "lncRNA\tno_lncRNA\tlncRNA_liftover\tno_lncRNA_liftover\tOne_to_one_orth\tone_to_half\tOne_to_many_orth\tMany_to_many_orth\tOne_to_none_orth\tOne_lncRNA_to_one_lncRNA\tOne_lncRNA_to_one_protein\tOne_lncRNA_to_one_other\n"
        lifted_spesi_data_lncRNA = lifted_spesi_data_lncRNA + "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(total_lncRNA_no, no_lncRNA2, lifted_lncRNA_no, not_lifted_lncRNA_no, lifted_one_one_lncRNA_no, lifted_one_half_lncRNA_no, lifted_one_many_lncRNA_no, lifted_many_many_lncRNA_no, lifted_one_none_lncRNA_no, lifted_one_one_lncRNA_lncRNA_no, lifted_one_one_lncRNA_pc_no, lifted_one_one_lncRNA_other_no)

        output_file_stat1 = os.path.join(statistic_2l, cleaned_filename2 + ".tsv")
        with open(output_file_stat1, 'w') as output_file_stat1:
            output_file_stat1.write(f"{lifted_spesi_data_lncRNA}\n") 
        subprocess.call("Rscript ./scripts/lncRNA_plot.R "+str(gene_level)+" > /dev/null 2>/dev/null", shell=True)

        lifted_spesi_data_pc = "pc\tno_pc\tpc_liftover\tno_pc_liftover\tOne_to_one_orth\tone_to_half\tOne_to_many_orth\tMany_to_many_orth\tOne_to_none_orth\tOne_pc_to_one_pc\tOne_pc_to_one_lncRNA\tOne_pc_to_one_other\n"
        no_pc2 = 0
        lifted_spesi_data_pc = lifted_spesi_data_pc + "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(total_pc_no, no_pc2, lifted_pc_no, not_lifted_pc_no, lifted_one_one_pc_no, lifted_one_half_pc_no, lifted_one_many_pc_no, lifted_many_many_pc_no, lifted_one_none_pc_no, lifted_one_one_pc_pc_no, lifted_one_one_pc_lncRNA_no, lifted_one_one_pc_other_no)

        output_file_stat2 = os.path.join(statistic_2p, cleaned_filename2 + ".tsv")
        with open(output_file_stat2, 'w') as output_file_stat2:
            output_file_stat2.write(f"{lifted_spesi_data_pc}\n") 
        subprocess.call("Rscript ./scripts/protein_plot.R "+str(gene_level)+" > /dev/null 2>/dev/null", shell=True)

###########################################################################################################################
###########################################################################################################################
###############################                                                        ####################################
###############################               Bi-Directional End                       ####################################                  
###############################                                                        ####################################
###########################################################################################################################
###########################################################################################################################
print("\ncompleted")
