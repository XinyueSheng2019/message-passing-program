# Percolate Program

## Code structure

1. percolate.h : The main header file.
  This file includes 5 structs: topology, neighbor, mapinfo, dynamicmaps and communication. Besides, it also contains many functions' declarations.
2. percolate.c : The main c file. This file includes the main function of this program and other essential functions:
```
  * divide_map(): For each process, after receving the whole map, it needs to extract the corresponding small map from the map.
  * send_smallmap(): send the small map to rank 0 process
  * combine_smallmap(): receive the smallmap from each processor and combine them into map(with extension)
  * drop_extension(): drop the extension (NL-L) part of the map (when L cannot be divided by the number of processes)
  * judge_percolate(): judge whether this map is percolated
```
3. initmap.c : This file only has one function init_map().
```
  * init_map(): This function can initialize the map, and set 0 to filled and expanded squares.
```
4. deliversmallmap.c : This file only has one function deliver_smallmap_to_old(). 
```
  * deliver_smallmap_to_old(): This function can deliver the smallmap in each process to a old map with width N+2 and length M+2, and then add halos to each old map.
```
5. percprocess.c : This file contains the loop of percolating each small map. 
```
  * percolate_process(): this is the major function of percolating. The following functions are called in this function.
  * swap_halos(): swap the old map' halos with its neighbors by using non-blocking communication methods.
  * extension_swap(): if the map is extended, some small maps need to swap the real halos.
  * compare_squares(): for each square, compare its neighbors' digits and refresh itself with the biggest digit.
  * get_sum_of_smallmap_digit(): calculate the sum of digits on smallmap.
```
6. arralloc.c : This file can be used to create a dynamic array.
7. uni.c : This file can generate random numbers.
8. percwritedynamic.c : This file can generate map.pgm according to the final map which is a dynamic array.


## Build the code on Cirrus
please load the mpt and intel-compiler-18:
```bash
module load intel-compiler-18 mpt
```
## Run the program
1. Makefile is included in this folder. Please input 'make' to compile the code:
```bash
make
```
2. After that, we could execute the code with n processes and a seed (0 ~ 900,000,000). You can randomly choose the number of processes and the seed. For example:
```bash
mpirun -n 16 ./percolate 1564
```
3. The PGM file named map can be generated. You could input this command to see it:
```bash
display map.pgm
```
4. If you want to delete all generated files, please input this command:
```bash
make clean
```



