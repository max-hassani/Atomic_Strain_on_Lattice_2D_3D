# Atomistic Strain mapping on a 2D/3D Lattice

This code can be used as a post-processing step to analyse the deformation in an atomistic simulation.

The code assumes that you have a series of dump files from MD simulation of LAMMPS. LAMMPS is an opensource parallel MD package. You can find more information on MD simulation and LAMMPS in their website (https://lammps.sandia.gov/).

**Building**  
This code includes four essential files: main.cpp, p_vec_ten.hh, p_vec_ten.cc and a Makefile. In a bash terminal, you can simply type `make` to build the code. Here, we assumed that you are apply this command in the directory where you have the four aformentioned files.

**Running the code**  
To run the code, first you need an input file. As an example, the main folder includes a file called "inputfile.txt". You need to pass to the code as an argument, the path to your input file.   
to run the code, you can simply insert the following command:  
` ./strain_map.X inputFile.txt
`   
For evaluating the atomistic strain, we need the positions of each atom. This is given in a dump file from lammps. Therefore, we need a series of dump files. Each dump file must have a name whic reflect its corresponding timestep, e.g., `dump.100000`.  
  
For each strain evaluation, one needs a pair of dump files. Therefore, at least we need two dump files. For our evaluation we need to give the initial timestep and the increment of timestep between dump files. For instance, our initial dump file is `dump.100000` and the next one is `dump.150000`. Therefore, our initial timestep would be 100000 and the increment is 50000. Each dump file must include the coordinate of the atoms (particles). The presumed format of the dump file is :

`atom_index atom_type xu yu zu`  


Here `xu yu zu` are the unwrapped position of the particles. The coordinates can also be normalized and replace by `xs ys zs`. For more information on this please refer to lammps documentation (https://lammps.sandia.gov/doc/dump.html). 

Before, the code starts its analysis, it needs several other parameters as well; the number of particles (atoms); The path to the dump files; information on the periodicity of your MD box; whether the coordinates of the particles are normalized by the box size or not; if user wishes a 2D mapping of the strain and displacement and in which cartesian plane; and finally whether you want the 3D data being saved or not.  

All these informations are included in the inputfile, path of which is given as an argument to the code.

**Structure of input file**  
The input file shoud include the following lines:  

`t_init=<value> 
t_step=<value> 
num_times=<value>  
NRPART=<value>
cut_off=<value>
dump_path=<value>
periodicity=<value>
normalized_coordinate=<value>
proj_plane=<value>
save_3D=<value>
`   

1) The value for t_init is an integer, which represents the initial timestep, you wish your analysis starts from. This will reflect itself in the name of the dump files, as mentioned above.  

2) The value for t_step is an integer, which inputs the code with the increment between two consecutive dump files.  

3) num_times need as value the number of evaluations of strains. This is also an integer, which is equal to the number of dump file minus one.  

4) NRPART represents the number of particles, therefore for its value an integer is expected.  

5) To evaluate the strain for each particle we need information about its neighboring particles. Here the value for the cut_off gives us the size of this neighborhood. It is recommend to choose this value as the first minimum in the radial distribution funtion, g(r). This value must be a float.   

6) dump_path, as its name suggests, needs as a value the path to the dump file. Here, I used a special format: `path_to_dump_file/name_*.dat`. The code replaces the astix with the timestep values from the dump file series. This value is a string.  

7) For periodicity, you need to pass a string, like `p s f`. Here we use LAMMPS convention for setting the periodicity of the box. Here `p` declares periodicity in x direction, while `s f` shows that in y and z direction we have a non-periodic boundary conditions.  

8) Normalized_coordinate expects as value, a `yes` or ` no`. Which inform the code, whether the coordinates are normalized or not.  

9) in the proj_plane, we set in which plane we need tht projection of the displacement/strain. Plausible values are `xy`, `xz`, `yz` or `none`. Here `none` declares that we do not desire any 2D projection. 


10) In save_3D, we declare whether it is necessary that the code save the 3D displacement/strain on the storage. The right value can be `yes` or `no`. The name of the saved files are set by the code itself. After the analysis you can move or rename the files.   


```python

```


```python

```


```python

```
