# HEP Forward and Backward Multiplicity


To run code: 
1. Change line 175 to the file path of your .txt containing list of root files (within the brackets)
'getListOfFiles("path/to/file.txt")'

2. Change line 465 to the file path of .root file to create
'TFile data_plot("/path/to/file/fb.root", "recreate");' 

3. In your OS's terminal, 'cd' into the folder containing the .txt and .cpp files and type 'root forward_backward.cpp' to compile and run code. 

4. In root session, type 'new TBrowser' and click in the GUI to view the plots generated in the root file. 
