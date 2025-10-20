# Usage: awk -v column={column_number} -v array={array} -f filter_results.awk input_file

# column_number is the column to filter
# array is a comma-separated string with one or more entries

# fields: chromosome, start, end, experiment_ID, TF_name, JASPAR_matrix_ID, score, strand

# field separator: \t

$(column) >= score {

  print

}

