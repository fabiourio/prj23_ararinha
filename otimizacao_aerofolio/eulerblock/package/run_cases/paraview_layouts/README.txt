The layout files will always open the solution.vtk file contained in the folder 02_single_run.

Use the following procedure to apply the state file:

- Open Paraview

- On the toolbar, select: File > Load State

- Select the paraview_state.py file

===========================================================

If you wish to open a different solution.vtk file, you should do the following:

- Open paraview_state.py

- Search for the line that contains "solutionvtk = LegacyVTKReader" and adjust the path to the desired solution.vtk file

- Open Paraview

- On the toolbar, select: File > Load State

- Select the paraview_state.py file

- A new window will appear. For the "Load State Data File Option" select "Choose File Names".

- Browse for a "solution.vtk" file in one of the run_case folders.

- Then click "OK".

- You should see the airfoil and CP curves.