import re
import pandas as pd

class database:
    def __init__(self, path = "./reference/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"):
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
    
    def parse_marker_file(self, path = "../../output/all/markers_btween_clusters.tsv", p_val_adj_threshold = 0.05):
        self.marker_dict = {}
        _ = pd.read_csv(path, sep = "\t")
        _ = _[_.p_val_adj <= p_val_adj_threshold]
        for cluster in set(_.cluster):
            self.marker_dict[cluster] = _[_.cluster == cluster].gene.to_list()
        print("--- Markers after filtering by p_val_adj({}):\n{}\n\n".format(p_val_adj_threshold, {k:len(v) for k,v in self.marker_dict.items()}))
    
    def filter_tf(self, gene):
        if not isinstance(gene, list):
            raise TypeError("gene parameter only takes list.")
        _ = [v for k,v in self.body.items() if k in gene]
        print(f"{len(_)} TFs found.")
        return [v for k,v in self.body.items() if k in gene]
    
    @staticmethod
    def parse_motif_name(motif_line):
        return motif_line.split()[1].split("_")[0]
    
    def wirte_meme(self, clusters = None, path = "./reference/"):
        if clusters is None:
            clusters = list(self.marker_dict.keys())
        
        if not isinstance(clusters, list):
            raise TypeError("clusters parameter only takes list.")
        
        for cluster in clusters:
            if cluster not in self.marker_dict.keys():
                print("Input cluster name", cluster, "is invalid.\n")
            print("--- Writing cluster", cluster, "to", path + str(cluster)+"_TF.meme")
            with open(path + str(cluster)+"_TF.meme", "w") as f:
                f.writelines(self.header)
                for tf_info in self.filter_tf(self.marker_dict[cluster]):
                    f.writelines(tf_info)
            print("Finished.\n")
    
    @staticmethod           
    def parse_uniProt_idmapping(path = "./idmapping_2024_04_17.tsv", gene = None):
        if not isinstance(gene, list):
            raise TypeError("gene parameter only takes list.")
        _ = pd.read_table(path)
        _ = _[["From", "Entry Name"]]
        # for simplicity, I'm returning all genes together without considering replacing meme TF name
        _.index = _.From
        _ = _.loc[list(set(gene) & set(_.index))]
        return list(set(gene + _.From.to_list() + _.loc[:,"Entry Name"].str.removesuffix("_HUMAN").to_list()))


if __name__=='__main__':
    d = database(path = "./HOCOMOCOv11_core_HUMAN_mono_meme_format.meme")
    activeGenes = ["Hendrikse2022_RLSVZ_activeGenes", "Hendrikse2022_RLVZ_activeGenes", "Hendrikse2022_RL_activeGenes"]

    for l in activeGenes:
        print(f"Processing {l}")
        with open(f"../activeGenes/{l}", "r") as g:
            genes = [line.strip() for line in g]

        with open(f"./filteredTFdb/{l}_TF.meme", "w") as f:
            f.writelines(d.header)
            for tf_info in d.filter_tf(gene = d.parse_uniProt_idmapping(path = "./idmapping_2024_04_17.tsv", gene = genes)):
                    f.writelines(tf_info)
            print("Finished.\n")










