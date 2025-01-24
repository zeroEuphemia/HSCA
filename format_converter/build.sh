cd src
make
cd ..
mv src/FormatConverter ./FormatConverter

cd CTW_Parser
g++ -O3 -std=c++11 parser.cpp -o ../ctw_parser
cd ..
