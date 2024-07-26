from dolfin import *
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

directory = "PlasticCorrectorResults"
if not os.path.exists(directory):
    os.makedirs(directory)
    print(f"Directory {directory} created")
else:
    print(f"Directory {directory} already exists")

def set_values_to_mesh():

    filename = f"PlasticCorrectorResults/p_plasticcorrector.pvd"
    file_p_plasticcorrector = File(filename)

    def MeshComposite():
        mesh = Mesh()
        filename = f"mesh_micro_10030.xdmf"
        print(f'loading mesh: {filename}')
        try:
            with XDMFFile(filename) as infile:
                infile.read(mesh)
        except:
            pass
        mvc = MeshValueCollection("size_t", mesh, 2) 
        mf = []
        return mesh, mf

    mesh, _ = MeshComposite()

    p_plasticcorrector_list = np.loadtxt(fname=f"./p_plasticcorrector.txt")
    num_time_steps = np.shape(p_plasticcorrector_list)[0]
    Z = FunctionSpace(mesh, "DG", 0)
    p_plasticcorrector = Function(Z)

    for i in range(num_time_steps):
        p_plasticcorrector.vector().set_local(p_plasticcorrector_list[i, :])
        p_plasticcorrector.rename("Plastic corrector p", "Plastic corrector p")
        file_p_plasticcorrector << (p_plasticcorrector, i)  # Write with time step
        print(f'Plastic corrector p set in mesh, time-step {i}')

set_values_to_mesh()


