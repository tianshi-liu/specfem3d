README
======

This folder contains several examples with increasing complexities.


To familiarize yourself with the package and CUBIT support,
we suggest to start with:

- **homogeneous_halfspace/**: a model with a single constant material description.

  Use the script `run_this_example.sh` in that directory to run the simulation,
  or execute the step-by-step instructions in the README from that directory.

More complex examples are provided in:

- **layered_halfspace/**: a model with different material layers used for benchmarking

- **waterlayered_halfspace/**: a model with fluid-solid coupling

- **tomographic_model/**: a model used together with of a tomographic file for the material description

- **noise_tomography/**: uses the homogeneous_halfspace model for noise kernel simulations

- **Mount_StHelens/**: a model with topography data included

- **homogeneous_poro/**: a model using poroelastic material


To familiarize yourself with the new in-house mesher xmeshfem3D (no CUBIT needed),
you can follow the examples in:

- **meshfem3d_examples/**


To familiarize yourself with simulations using CPML absorbing boundaries,
you can follow the examples in:

- **CPML_examples/**


To familiarize yourself with simulations using dynamic/kinematic rupture models,
you can follow the examples in:

- **fault_examples/**


Please **consider submitting your own example** to this package!



##### ADDITIONAL NOTES:
- Because of the poroelastic implementation, even in an acoustic/elastic case,
    the file nummaterial_poroelastic_file must be present in the MESH/ directory.
    Need to had a dummy file in cubit2specfem3d.py [CM - Oct. 2011]

    For now it can be copied from `homogeneous_poroelastic/MESH_coarse`.

- For the examples using the CUBIT_GEOCUBIT python scripts, please see the `CUBIT_GEOCUBIT/README.md` file
    about how to set up the CUBIT (now called TRELIS) and GEOCUBIT meshing packages.


