cd generate
sh build.sh
cd ..

cd optimize
sh build.sh
cd ..

cd format_converter
sh build.sh
cd ..

mv ./generate/SamplingCA/SamplingCA ./SamplingCA
mv ./generate/Generator ./Generator
mv ./optimize/Optimizer ./Optimizer
