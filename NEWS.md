# nemoRsem 0.2.0.0 
<span style="color:grey;">10/04/2024</span>

* `scheme_interp_slice` now returns a scheme to extract an array layer if the target depth already exists in the dataset.
* `get_*_slabr` functions now have a collapse_days argument which lets you average across timesteps within a file. 
* `process_array` is now used within `get_*_slabr` functions for exception handling of array dimensions. 
* `NE_volume_summary` no longer uses hard-coded column names.

# nemoRsem 0.1.0.0 
<span style="color:grey;">16/05/2023</span>

* Copied nemomedusR package and began updating variables and documentation for the NEMO-ERSEM model output
