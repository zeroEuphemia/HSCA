echo "c now call SamplingCA to generate an initial 2-wise CA"
./SamplingCA/SamplingCA -seed $1 -input_cnf_path $2 -output_testcase_path ./SamplingCA/$3
echo "c now call ScalableCA to generate a 3-wise CA"
./ScalableCA -seed $1 -input_cnf_path $2 -init_CA_file_path ./SamplingCA/$3 -output_testcase_path $3 -strength $4 -group_file_path $5 -L 5000 -use_group 1 -opt_method 1 -use_weight 1
rm ./SamplingCA/$3
