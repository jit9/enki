# ptsrc-cat
Point source catalog code
Written by Toby Marriage, Megan Gralla... (used by Rahul, Kevin, Heather)

Before running this code, ensure that you have flipper installed. 
If you want to use Enki maps then you may want to have ndflip too.

To run this code you will need intensity and weight maps (from Enki or Ninkasi) 
and a beam datafile


1. Clone the point source code from github:
git clone https://github.com/ACTCollaboration/ptsrc-cat.git


2. Edit your bashrc or bash_profile to point to the point source code so you
can call it from elsewhere
(there are neater ways of doing this, should we implement them?)
In the following, POINT_SOURCE_DIR is the path to the directory containing the
bin and python folders.

  export POINT_SOURCE_DIR=path to where point source code is
  export PATH=$PATH:$POINT_SOURCE_DIR/bin
  export PYTHONPATH=$PYTHONPATH:$POINT_SOURCE_DIR/python
  export PYTHONPATH=$PYTHONPATH:$POINT_SOURCE_DIR/bin


3. Create your input maps. To make a catalog the code needs the total
intensity map in Jy/sr and a weight map. Depending on the input maps used
there are different ways to do this.

If you are using maps that do not have the data split then for Enki, use ndflip
to separate out the intensity map and weights, convert the intensity to Jy/sr
and save. An example of this is given in make_input_maps_enki.py. For Ninkasi
maps you may need to add together the I, beam_test and input_used maps to get
the total intensity.

You may also need to add together maps from a four way split. I (Heather) have
not done this but the approach outlined to me by Kevin was

  i. add together the CMB I and beam_test maps using eg add_input_beam_test_and_100.py

  ii. Run weightedCoadd which uses weightedCoadd.dict

  iii. Run get_maps - not sure if this does anything important or just renames files

  iv. Run convertToJyPerSr to change from micro Kelvin to Janskys per steradian

Whatever method you use, you need an input intensity map that is in Jy per sr
and a weight map to use as input to makeCatalogMaster by the end of this step.


4. Modify the input dictionaries:
makeCatalogMaster.dict
  - ensure that map and weight point to the input intensity and weight
    maps that you have created
  - check that the region(s) of the map where you want to look for sources
    are defined in fielddefs
  - specify the region(s) you are interested in in fieldnames
  - ensure that lowestWeightAllowed is a reasonable cut given your input
    weight map (the number may be higher if your weight map is just a
    hit count and lower if detector noise is taken into account)
  - select your initial and final signal to noise ratio thresholds

makeCatalog.dict
  - set signalTransform1D to point to the beam datafile (in beams)
  - set extraFilter1D to point to an optional additional (eg high pass) filter
  - if you want fits files for each source then set writeSubmaps=True
    (caution: if you have lots of sources these files will take up a lot
    of space and may slow your code down)

makeTemplateFromCatalog.dict
  - set templateMap to ../data_YOUR_FIELD_NAME.fits
    where YOUR_FIELD_NAME is the field name from makeCatalogMaster.dict
  - set file in templates to point to the beam datafile in beams


5. Create the catalog by running:
makeCatalogMaster makeCatalogMaster.dict

  - if the catalog was created successfully the code will produce a folder
    with the name you specified in fieldnames in makeCatalogMaster.dict.
    Inside this folder is init, which has the initial catalog of high signal to
    noise sources, and main which has the rest of the sources above the signal
    to noise threshold.

Note: do not despair if your code crashes with the error
IOError: [Errno 2] No such file or directory: 'FIELDNAME/main/catalogs/catalog.pickle.dict'
To find the actual error, navigate to FIELDNAME/init and check makeCatalog.err
or makeTemplate.err to see where and why things crashed.
