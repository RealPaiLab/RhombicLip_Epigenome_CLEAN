### get gene expressed in > 1% cells in Hendrikse RL VZ/SVZ 
Rscript --no-save ../Hendrikse_2022.R

### remove lowly expressed genes from TF database
python3 filter_tf.py
