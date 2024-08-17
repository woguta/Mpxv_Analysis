# Mpxv_Anlaysis
This an analaysi of Mpox virus sequences

```
#!/bin/bash

# Input and output file paths
input_file="./mpox_august/mpox_433_all.fasta"
output_file="./mpox_august/mpox_433_all.clean.fasta"

# Remove sequence length from headers and write to output file
sed 's/\/[0-9]\+-[0-9]\+//' "$input_file" > "${output_file}.tmp1"

# Replace special characters in headers with underscores and ensure proper FASTA format
awk '
    BEGIN {header_num=1}
    /^>/ {
        header=$0
        gsub(/[|\/\-]/, "_", header)
        print header
    }
    !/^>/ {print}
' "${output_file}.tmp1" > "${output_file}.tmp2"

mv "${output_file}.tmp2" "$output_file"
rm "${output_file}.tmp1"

# Ensure the final fasta file has a proper FASTA header
if ! grep -q "^>" "$output_file"; then
    echo "ERROR!!! The file ${output_file} does not have a proper FASTA header. Adding default headers."
    awk '
        BEGIN {header_num=1}
        !/^>/ {
            print ">sequence_" header_num
            header_num++
        }
        {print}
    ' "$output_file" > "${output_file}.tmp"
    mv "${output_file}.tmp" "$output_file"
    echo "Default headers added to ${output_file}"
else
    echo "The file ${output_file} has been cleaned and is in proper FASTA format"
fi
```

