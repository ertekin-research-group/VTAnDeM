VTAnDeM
=======

VTAnDeM (Visualization Toolkit for Analyzing Defects in Materials) is a post-processing plotting toolkit for DFT calculations 
of defects in semiconductors. For given materials, the script imports DFT data provided by users, and produces phase stability 
diagram, defect diagrams at different chemical environments, and carrier concentrations. 


Getting Started
---------------
- Use Compound_Data_Extract.py to generate CSV file needed by the tool
- Choose two elements (First and Second Species) as the axis of phase diagram
- Press ‘Generate Plot!’ Button, use slide bar in the bottom to change chemical    environment
- Pick up a point in the shaded stability region
- Presse ‘Defect Diagram’ Button
- Choose different point in the stability region to explore defect diagrams


User-provided data input:
-------------------------
DFT calculation results (POSCARs and OUTCARs) for: 
-	bulk chemical constituents;
-	competing compounds;
-	main compound.
generated by:
-	scf calculation for bulk and defect structures.

Powered by
PyQT5
