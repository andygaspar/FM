rm fast_me_

gcc  -c test_main.c main.c distance.c bme.c bNNI.c gme.c graph.c heap.c inputs.c interface_options.c interface_utilities.c newick.c NNI.c p_utils.c SPR.c traverse.c random.c  utils.c -fopenmp
gcc -w -o fast_me_  distance.o bme.o bNNI.o gme.o graph.o heap.o inputs.o interface_options.o interface_utilities.o newick.o NNI.o p_utils.o SPR.o traverse.o utils.o random.o  main.o test_main.o -lm -fopenmp -lpthread

rm *o

cd ..
./src/fast_me_ -i mat.mat -m b -n -s