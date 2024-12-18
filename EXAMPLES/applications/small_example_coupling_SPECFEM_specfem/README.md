# Couple with injection - example of hybrid method run Specfem/Specfem

This is a very small example to test the coupling between SPECFEM3D and itself.
In other words, it runs a coarse simulation to create an initial wavefield for a distant event,
and then re-injects that wavefield at the coupling boundaries of a small-scale local simulation.

This coupling example uses the following steps:
1. create the coupling boundary points of the small-scale local mesh by running the mesher
   with the small-scale local setup in folder `DATA.local/`

2. run the coarse simulation using DATA.coarse/ to store the wavefield at the boundary points
   into a separate folder `COUPLING_FILES/`

3. run the small-scale local simulation with injecting the wavefield at the boundary points

The default setup uses 4 MPI processes.

To launch the full coupling example with the above 3 steps, type:
```
> run_this_example.sh
```

notes:

* The coupling example here between coarse and small-scale mesh uses a flat top surface.
  This allows the coarse simulation to locate all coupling points at the exact locations.
  The resulting seismograms for stations inside the small-scale region should therefore become (almost) identical.

  You can compare the seismograms in the output folder `OUTPUT_FILES/` from the small-scale local simulation, with the
  seismograms in `OUTPUT_FILES.coarse/` from the coarse simulation.

* In case coarse and small-scale simulation would use realistic topographic interfaces, there might be a difference
  in topography resolution and thus coupling points from the small-scale coupling boundary could be located slightly
  outside the coarser mesh. Locating those coupling points would move them closer to the nearest coarse mesh elements.
  Storing the wavefield and re-injecting it would then lead to smaller artefacts.

  To work-around this, one would have to make sure that the top of the small-scale mesh surfaces match exactly
  the top of the coarse mesh surfaces at the injection boundaries.
  This would either require to tweak the topo interfaces for the in-house mesher or use an external mesher to resolve this.
  For "smoother" incoming wavefields and small discrepencies between surface points, the injected wavefield however
  might still be good enough.

* To save storage space, we use `NTSTEP_BETWEEN_OUTPUT_SAMPLE == 2` in the `Par_file` of the coarse simulation.
  This seems to provide enough time samplings for the injection. One would have to test this out a bit for other examples.
  In general, the storage of the wavefield components (velocity & traction) over all coupling boundary points and over all time steps
  can lead to significant file sizes.

  In any case, the stored wavefield will be (time) interpolated during the preparation stage of the small-scale local simulation,
  as the time steps `DT` of the coarse and small-scale simulations are likely different. Time interpolation is done by a
  Catmull-Rom scheme that provides cubic interpolation between time samples. This scheme seems to exhibit a good balance between accuracy
  and speed for these kind of simulations. Thus, storing the wavefield at lower sampling rates can be fine for many applications.



