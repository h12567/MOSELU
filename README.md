# MOSELU

MOlecule Structure ELucidation (MOSELU)

A benchmarking platform for reconstructing the molecule structure given observed mass spectrum.

## Reading list

1. https://www.cpp.edu/~psbeauchamp/pdf_book/MS_chapter.pdf
2. https://www.ics.uci.edu/~dock/manuals/DaylightTheoryManual/theory.smarts.html



## Requirements

https://www.rdkit.org/docs/Install.html
highly recommend creating a conda environment as they have suggested
$ conda create -c rdkit -n my-rdkit-env rdkit

pip freeze for your convenience, no need stick too closely, install specific version only when there is exception.

antlr4-python3-runtime==4.9.1
certifi==2020.12.5
importlib-metadata==3.4.0
jsonpickle==1.5.1
mkl-fft==1.2.0
mkl-random==1.1.1
mkl-service==2.3.0
numpy @ file:///tmp/build/80754af9/numpy_and_numpy_base_1603479632437/work
olefile==0.46
pandas==1.2.1
Pillow @ file:///tmp/build/80754af9/pillow_1609786740831/work
python-dateutil @ file:///home/ktietz/src/ci/python-dateutil_1611928101742/work
python-igraph==0.8.3
pytz @ file:///tmp/build/80754af9/pytz_1612215392582/work
scipy==1.6.0
six @ file:///tmp/build/80754af9/six_1605205313296/work
texttable==1.6.3
tqdm==4.56.0
typing-extensions==3.7.4.3
zipp==3.4.0



## Generating subset of dataset

$ cd MOSELU
$ python nist_db_helpers/prepare_train_dataset.py

run it like this, some directories are hardcoded
"list assignment index out of range" errors can be ignored
Expect the output files to be of length 4882, assuming no parameters was changed.

change line 171-175 of prepare_train_dataset.py for your needs.
func_groups = ["ester", "ether"]
allow_molecules = ["C", "H", "O", "N"]
max_constraint = [("C", 11)]
possible_bonds = [BondType.SINGLE, BondType.DOUBLE, BondType.TRIPLE]
max_atoms = 13

For allowed func_groups, check out ./data/smarts.json, don't edit this file! To understand it, check out reading #2.
max_constraint is basically the maximum number of ? atom we allow in the molecule. In this case we only allow 11 carbons at most.
possible_bonds should be straight forward, if I'm not wrong can have quadruple, and aromatic....
max_atoms DOES not include hydrogen atoms in the count. So it is possible to have a C11H24.

note that original dataset has no atom-ordering,
so order is implemented in line 117's extract_adj_matrix_and_order_vertices().
basically what happens is for a given ester,

A
((4, 2, 3, 1, 0),)
['C', 'C', 'O', 'O', 'C']
['C', 'O', 'C', 'O', 'C', 'C', 'O', 'N', 'C', 'C', 'O']
B
[[0. 1. 0. 0. 0.]
 [1. 0. 2. 1. 0.]
 [0. 2. 0. 0. 0.]
 [0. 1. 0. 0. 1.]
 [0. 0. 0. 1. 0.]]
[0 0 2 2 0]

Original unordered is line 66 of this file,
Ester functional group is detected, with fixed order C C O O C (line 65) with fixed bond pattern (line 68-72)
Line 64 is the index of the stuff in line 65 in the original unordered molecule in line 66.

For the remaining atoms which are not part of the functional group, they will be arranged in increasing index.
Read MassPreserveSmilesOrdering in ordering.py for exact implementation
So ordered output will be [0 0 2 2 0 0 0 0 2 2 3] or CCOOCCCCOON



## Data files explanation

4 files produced

###1

>>> np.load('func_groups.npy',allow_pickle=True).shape
(4882,)
>>> np.load('func_groups.npy',allow_pickle=True)[0]
'ester'

functional group label

###2

>>> np.load('vertex_arr_sort_svd.npy',allow_pickle=True).shape
(4882,)
>>> np.load('vertex_arr_sort_svd.npy',allow_pickle=True)[0]
array([0, 0, 2, 2, 0, 0, 0, 0, 0])

0 means carbon, 2 means oxygen. number is determined by indexes of allow_molecules in prepare_train_dataset.py

###3

>>> np.load('mol_adj_arr_sort_svd.npy', allow_pickle = True).shape
(4882, 13, 13)
>>> np.load('mol_adj_arr_sort_svd.npy', allow_pickle = True)[0]
array([[0., 1., 0., 0., 0., 2., 1., 0., 0., 0., 0., 0., 0.],
       [1., 0., 2., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [1., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])

Conditioned on the order of the vertex_arr_sort_svd.npy file, this one shows the single/double/triple bond between some pair of atoms. Sufficiently indicative of structure.

###4

>>> np.load('msp_arr.npy', allow_pickle = True).shape
(4882, 800)
>>> np.load('msp_arr.npy', allow_pickle = True)[0]
array([  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   9.,  39., 592.,  11.,  20.,  72.,   0.,   0.,   0.,
         0.,   0.,   0.,   5.,  81., 448.,  98., 133.,  24.,  30.,   2.,
         6.,   0.,   0.,   0.,  14.,  39., 394.,  70., 999.,  96., 125.,
        33.,  62.,   0.,   0.,   0.,   5.,  46.,  82.,  40., 323.,  94.,
       160.,  22.,  19.,   6., 349.,  12.,   4.,   5.,   9.,   3.,  21.,
        18., 174., 187., 715.,  46.,  17.,   6., 242.,  27.,   0.,   0.,
         0.,   0.,   3.,   0., 140.,  17.,  13.,   0.,  50.,  19.,  34.,
         4.,   0.,   0.,   0.,   0.,   0.,   0.,  39., 230., 307.,  23.,
        10.,   5.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0., 142.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   6., 229.,  19.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.])

Each mass-spec is a vector of size 800, with mass/charge of 0 to 799. Max peak is by definition 999, everything else is relative to the max peak. So this molecule has a max peak at mass/charge of 41.