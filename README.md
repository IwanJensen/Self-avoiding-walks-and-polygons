# Self-avoiding walks and polygons on the square lattices
Series expansions for SAW and SAP problems on the square lattice <pre>
sqsaw.ser          Number of SAW 
sqsaw_EE.ser       End-to-end distance of SAW
sqsaw_RG.ser       Radius of gyration of SAW
sqsaw_MD.ser       Mean monomer distance of SAW
sqsaw_end_distibution.ser   Complete end-to-end distribution of SAWs of length n starting at (0,0) and ending at (x,y). Format of data is x, y, n, c_n(x,y). By symmetry we only list data for x>=1, y>=0, y<=x.
sqsaw_wormss.ser   Number pf worms (y=0 data from sqsaw_end_distibution.ser).

sqsap_perim.ser    Number of SAP by perimeter
sqsap_area.ser     Number of SAP by area
sqsap_RG.ser       Radius of gyration of SAP
