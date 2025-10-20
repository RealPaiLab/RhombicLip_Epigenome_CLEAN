import os
import re
import pandas as pd

class database:
    """
    Parse HOCOMOCOv12 meme format database and perform basic motif filtering.
    """
    def __init__(self, path = "./reference/H12CORE_meme_format.meme"):
        self.__read_database(path)
        print(f"--- Read in {len(self.body)} TFs in total.\n\n")
    
    def __read_database(self, path):
        with open(path, "r") as f:
            _ = f.readline()
            header = [_]
            while not re.match("A.+C.+G.+T.+", _):
                _ = f.readline()
                header.append(_)
            header.append(f.readline())
            
            body = {}
            for line in f:
                if line.startswith("MOTIF"):
                    key = self.parse_motif_name(line)
                    value = [line]
                else:
                    while line != "\n":
                        value.append(line)
                        line = f.readline()
                    body[key] = value + ["\n"]
            self.header = header
            self.body = body
    
    def filter_tf(self, motif_names):
        if not isinstance(motif_names, list):
            raise TypeError("gene parameter only takes list.")
        _ = [v for k,v in self.body.items() if k in motif_names]
        print(f"{len(_)} TFs found.")
        return _
    
    @staticmethod
    def parse_motif_name(motif_line):
        return motif_line.split()[1]
    
    @staticmethod
    def parse_annoation_jsonl(path = "./reference/H12CORE_annotation.jsonl", genes = None):
        """ 
        Parse HOCOMOCO annotation jsonl and return a list of motif names of the provided genes
        """
        if not isinstance(genes, list):
            raise TypeError("genes parameter only takes list.")
            
        _ = pd.read_json("./reference/H12CORE_annotation.jsonl", lines = True)
        _ = pd.json_normalize(_.to_dict(orient='records')) # to parse nested structure
        
        # 1. Be aware of that original_motif.tf is not always the same as masterlist_info.tf
        # masterlist_info.tf is the one should be used
        # 2. If there are genes whose symbol and synonyms do not agree in expression, some genes 
        # may need to be removed will be kept. e.g. IRF4 and MUM1 are synonyms. IRF4 is lowly expressed
        # in RL, but MUM1 isn't. Then any TF related to MUM1 will be kept. This is not very common and
        # keeping the gene can avoid false negative. 
        rule = _.apply(lambda x: \
               len(set([x["masterlist_info.species.HUMAN.gene_symbol"]] + \
                       x["masterlist_info.species.HUMAN.gene_synonyms"]) & \
                   set(genes)) > 0, 
               1)
        return _.loc[rule, "name"].unique().tolist()



if __name__=='__main__':
    d = database(path = "/.mounts/labs/pailab/src/transcription-factors/HOCOMOCOv12/H12CORE_meme_format.meme")
    activeGenes = ["Hendrikse2022_RLSVZ_activeGenes", "Hendrikse2022_RLVZ_activeGenes", "Hendrikse2022_RL_activeGenes"]

    filteredTFdb_dir = "/data/xsun/db/meme/HOCOMOCOv12/filteredTFdb"
    if not os.path.exists(filteredTFdb_dir):
        os.makedirs(filteredTFdb_dir)
        print(f"--- {filteredTFdb_dir} created to store filtered meme db.\n")
    else:
        print(f"--- {filteredTFdb_dir} exists. The output will be written under this directory.\n")

    for l in activeGenes:
        print(f"--- Processing {l}")
        with open(f"/data/xsun/db/meme/HOCOMOCOv12/activeGenes/{l}", "r") as g:
            genes = [line.strip() for line in g]
        
        with open(f"{filteredTFdb_dir}/{l}_TF.meme", "w") as f:
            f.writelines(d.header)
            for tf_info in d.filter_tf(motif_names = database.parse_annoation_jsonl(path = "/.mounts/labs/pailab/src/transcription-factors/HOCOMOCOv12/H12CORE_annotation.jsonl", genes = genes)):
                    f.writelines(tf_info)
            print("--- Finished.\n")




