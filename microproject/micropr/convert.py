from pyvista import read
grid = read("CC=CCC=CCC=CCCCCCCCCC(=O)O_mo_0.cube")
grid.save("output_lenol.vti")