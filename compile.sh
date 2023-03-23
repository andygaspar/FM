

cd src
rm fast_me_

g++  -c test_main.cc fastme.c initialiser.c distance.c bme.c bNNI.c gme.c graph.c heap.c inputs.c interface_options.c interface_utilities.c newick.c NNI.c p_utils.c SPR.c traverse.c random.c  utils.c -fopenmp
g++ -w -o fast_me_  distance.o bme.o bNNI.o gme.o graph.o heap.o inputs.o interface_options.o interface_utilities.o newick.o NNI.o p_utils.o SPR.o traverse.o utils.o random.o initialiser.o fastme.o test_main.o -lm -fopenmp -lpthread

rm *o
echo "$PWD"



cd ..
./src/fast_me_ -i mat.mat -m b -n -s -f 17
./fastme -i mat.mat -m b -n -s -u init_topology.nwk -f 17

